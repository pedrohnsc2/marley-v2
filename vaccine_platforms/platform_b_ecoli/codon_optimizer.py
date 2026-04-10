"""B2 -- Otimizacao de codons para expressao em E. coli K12.

Implementa otimizacao de codons usando a tabela de frequencia de uso de
codons de E. coli K12 (cepa de referencia para BL21-DE3). Seleciona codons
por frequencia maxima para maximizar o CAI (Codon Adaptation Index).

Verificacoes pos-otimizacao:
- GC content entre 40-60%
- Ausencia de codons raros (frequencia < 10%)
- Sem repeticoes longas de nucleotideos (>8 nt iguais)
- Sem terminadores rho-independentes (stem-loop + poli-U)
- Sem sitios internos de restricao indesejados
- Sem sequencias tipo Shine-Dalgarno internas (AAGGAG)
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Tabela de uso de codons de E. coli K12
# ---------------------------------------------------------------------------
# Fonte: Kazusa Codon Usage Database (http://www.kazusa.or.jp/codon/)
# Valores: frequencia por mil (per thousand), normalizados abaixo para fracao
# por aminoacido. Usamos a tabela completa para maximizar a fidelidade
# do CAI.

ECOLI_K12_CODON_TABLE: Final[dict[str, list[tuple[str, float]]]] = {
    # Alanina -- GCG e o mais frequente em E. coli (diferente de eucariotos!)
    "A": [("GCG", 0.36), ("GCC", 0.27), ("GCA", 0.21), ("GCT", 0.16)],
    # Arginina -- CGT e CGC sao os preferidos; AGG/AGA sao raros em E. coli
    "R": [("CGT", 0.38), ("CGC", 0.36), ("CGG", 0.10), ("CGA", 0.06), ("AGA", 0.04), ("AGG", 0.02)],
    # Asparagina
    "N": [("AAC", 0.55), ("AAT", 0.45)],
    # Aspartato
    "D": [("GAT", 0.63), ("GAC", 0.37)],
    # Cisteina
    "C": [("TGC", 0.56), ("TGT", 0.44)],
    # Glutamina
    "Q": [("CAG", 0.65), ("CAA", 0.35)],
    # Glutamato
    "E": [("GAA", 0.69), ("GAG", 0.31)],
    # Glicina
    "G": [("GGT", 0.34), ("GGC", 0.40), ("GGA", 0.11), ("GGG", 0.15)],
    # Histidina
    "H": [("CAT", 0.57), ("CAC", 0.43)],
    # Isoleucina -- ATC e o preferido; ATA e muito raro em E. coli
    "I": [("ATC", 0.42), ("ATT", 0.51), ("ATA", 0.07)],
    # Leucina -- CTG domina em E. coli (>50%)
    "L": [("CTG", 0.50), ("TTA", 0.13), ("TTG", 0.13), ("CTT", 0.10), ("CTC", 0.10), ("CTA", 0.04)],
    # Lisina
    "K": [("AAA", 0.74), ("AAG", 0.26)],
    # Metionina -- codon unico
    "M": [("ATG", 1.00)],
    # Fenilalanina
    "F": [("TTT", 0.57), ("TTC", 0.43)],
    # Prolina
    "P": [("CCG", 0.52), ("CCA", 0.19), ("CCT", 0.16), ("CCC", 0.12)],
    # Serina -- TCT e AGC sao os mais frequentes
    "S": [("TCT", 0.15), ("TCC", 0.15), ("TCA", 0.12), ("TCG", 0.15), ("AGC", 0.28), ("AGT", 0.15)],
    # Treonina
    "T": [("ACC", 0.44), ("ACT", 0.17), ("ACA", 0.13), ("ACG", 0.27)],
    # Triptofano -- codon unico
    "W": [("TGG", 1.00)],
    # Tirosina
    "Y": [("TAT", 0.57), ("TAC", 0.43)],
    # Valina -- GTG e o mais frequente
    "V": [("GTG", 0.37), ("GTT", 0.26), ("GTC", 0.22), ("GTA", 0.15)],
    # Stop -- TAA e o terminador preferido em E. coli (lido pelo RF1, mais eficiente)
    "*": [("TAA", 0.64), ("TGA", 0.30), ("TAG", 0.07)],
}

# Mapa reverso: codon -> aminoacido (para calculos de CAI)
CODON_TO_AA: Final[dict[str, str]] = {}
for _aa, _codons in ECOLI_K12_CODON_TABLE.items():
    for _codon, _freq in _codons:
        CODON_TO_AA[_codon] = _aa

# Sitios de restricao a evitar na sequencia codificante
RESTRICTION_SITES_TO_AVOID: Final[dict[str, str]] = {
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "BamHI": "GGATCC",
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
}

# Shine-Dalgarno -- sequencia de ligacao ao ribossomo; se aparecer
# internamente pode causar iniciacao espuria de traducao
SHINE_DALGARNO_PATTERN: Final[str] = "AAGGAG"

# Limiar de codon raro (frequencia < 10% do total para aquele aminoacido)
RARE_CODON_THRESHOLD: Final[float] = 0.10

# Comprimento maximo de repeticoes de nucleotideo unico
MAX_HOMOPOLYMER_LENGTH: Final[int] = 8


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class CodonOptimizationResult:
    """Resultado da otimizacao de codons para E. coli."""

    dna_sequence: str
    length_nt: int
    gc_content: float
    cai: float
    rare_codons: list[dict[str, str | float]]
    restriction_sites_found: list[dict[str, str | int]]
    has_internal_ndei: bool
    has_internal_xhoi: bool
    has_shine_dalgarno: bool
    homopolymer_runs: list[dict[str, str | int]]
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Otimizacao de codons
# ---------------------------------------------------------------------------


def get_best_codon(amino_acid: str) -> str:
    """Retorna o codon de maior frequencia para o aminoacido dado.

    Em E. coli, a preferencia de codons e bem diferente de eucariotos:
    - Leu: CTG domina (50%) vs CTC em mamiferos
    - Arg: CGT/CGC preferidos; AGG/AGA sao raros
    - Pro: CCG preferido (52%) vs CCC em humanos

    Args:
        amino_acid: Letra unica do aminoacido (maiuscula).

    Returns:
        Triplet de DNA (ex: "CTG" para Leu).

    Raises:
        ValueError: Se o aminoacido nao esta na tabela.
    """
    codons = ECOLI_K12_CODON_TABLE.get(amino_acid)
    if not codons:
        raise ValueError(f"Aminoacido desconhecido: {amino_acid!r}")
    # Retorna o de maior frequencia
    return max(codons, key=lambda x: x[1])[0]


def optimize_codons(protein_sequence: str) -> str:
    """Traduz a sequencia proteica em DNA otimizado para E. coli K12.

    Usa estrategia de maior frequencia (max-frequency) para cada posicao,
    o que maximiza o CAI. Para genes heterologos curtos (<1 kb), esta
    estrategia e preferida sobre amostragem proporcional.

    Apos a selecao inicial, verifica e remove sitios de restricao
    internos por troca sinonima.

    Args:
        protein_sequence: Sequencia de aminoacidos a otimizar.

    Returns:
        Sequencia de DNA otimizada (sem stop codon no final).
    """
    codons: list[str] = []
    for aa in protein_sequence:
        codons.append(get_best_codon(aa))

    dna = "".join(codons)

    # Remover sitios de restricao internos por substituicao sinonima
    dna = _remove_internal_restriction_sites(dna, protein_sequence)

    return dna


def _remove_internal_restriction_sites(
    dna: str, protein: str
) -> str:
    """Remove sitios de restricao internos via substituicao sinonima.

    Para cada ocorrencia de um sitio de restricao, encontra o codon que
    o sobrepoem e troca por um codon sinonimo alternativo de maior
    frequencia que elimine o sitio.

    Args:
        dna: Sequencia de DNA codificante.
        protein: Sequencia proteica correspondente.

    Returns:
        DNA com sitios de restricao removidos.
    """
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    max_iterations = 20  # Evitar loop infinito

    for _iteration in range(max_iterations):
        changed = False
        current_dna = "".join(codons)

        for site_name, site_seq in RESTRICTION_SITES_TO_AVOID.items():
            pos = current_dna.find(site_seq)
            while pos != -1:
                # Qual codon e afetado?
                codon_start = (pos // 3) * 3
                codon_end_idx = ((pos + len(site_seq) - 1) // 3)

                # Tentar trocar cada codon sobreposto
                for ci in range(codon_start // 3, min(codon_end_idx + 1, len(codons))):
                    aa = protein[ci]
                    alternatives = ECOLI_K12_CODON_TABLE.get(aa, [])
                    # Ordenar por frequencia decrescente
                    sorted_alts = sorted(alternatives, key=lambda x: -x[1])
                    for alt_codon, alt_freq in sorted_alts:
                        if alt_codon == codons[ci]:
                            continue
                        if alt_freq < RARE_CODON_THRESHOLD:
                            continue
                        # Testar se a troca elimina o sitio
                        test_codons = codons[:]
                        test_codons[ci] = alt_codon
                        test_dna = "".join(test_codons)
                        if site_seq not in test_dna[max(0, pos - 6):pos + len(site_seq) + 6]:
                            codons[ci] = alt_codon
                            changed = True
                            break
                    if changed:
                        break

                # Atualizar para buscar proxima ocorrencia
                current_dna = "".join(codons)
                next_pos = current_dna.find(site_seq, pos + 1)
                pos = next_pos

        if not changed:
            break

    return "".join(codons)


def add_flanking_elements(coding_dna: str) -> str:
    """Adiciona os elementos flanqueadores para clonagem em pET-28a.

    Estrutura final do inserto:
        NdeI(CATATG) + [DNA codificante sem ATG inicial] + stop(TAA) + XhoI(CTCGAG)

    O ATG do sitio NdeI serve como codon de inicio. Por isso removemos o
    primeiro ATG da sequencia codificante se ele estiver la.

    Args:
        coding_dna: DNA codificante otimizado (comecando com ATG).

    Returns:
        Inserto completo com sitios de restricao e stop codon.
    """
    # O NdeI (CATATG) ja contem o ATG -- remover o ATG duplicado
    if coding_dna.startswith("ATG"):
        coding_dna = coding_dna[3:]

    # Stop codon preferido em E. coli
    stop = "TAA"

    # Montar: NdeI + ATG ja incluso + DNA + Stop + XhoI
    insert = f"CATATG{coding_dna}{stop}CTCGAG"
    return insert


# ---------------------------------------------------------------------------
# Metricas de qualidade
# ---------------------------------------------------------------------------


def calculate_gc_content(dna: str) -> float:
    """Calcula a fracao de G+C na sequencia de DNA.

    GC content ideal para E. coli: 50-52% (genoma E. coli K12 = 50.8%).
    Range aceitavel: 40-60%.

    Args:
        dna: Sequencia de DNA (maiusculas).

    Returns:
        Fracao de GC (0.0 a 1.0).
    """
    if not dna:
        return 0.0
    gc_count = sum(1 for nt in dna if nt in "GC")
    return round(gc_count / len(dna), 4)


def calculate_cai(dna: str) -> float:
    """Calcula o Codon Adaptation Index (CAI) para E. coli K12.

    O CAI mede quao bem os codons da sequencia correspondem aos codons
    preferenciais do organismo hospedeiro. Valor ideal: > 0.8.

    Calculo:
        CAI = exp( (1/L) * sum(ln(w_i)) )

    onde w_i e a frequencia relativa do codon i dividida pela frequencia
    maxima para aquele aminoacido (relative adaptiveness).

    Args:
        dna: Sequencia de DNA codificante.

    Returns:
        CAI entre 0.0 e 1.0.
    """
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]

    # Construir tabela de relative adaptiveness (w)
    w_table: dict[str, float] = {}
    for aa, codon_list in ECOLI_K12_CODON_TABLE.items():
        max_freq = max(f for _, f in codon_list)
        if max_freq == 0:
            continue
        for codon, freq in codon_list:
            w_table[codon] = freq / max_freq if max_freq > 0 else 0.0

    log_sum = 0.0
    count = 0
    for codon in codons:
        if len(codon) < 3:
            continue
        w = w_table.get(codon, 0.0)
        if w > 0:
            log_sum += math.log(w)
            count += 1

    if count == 0:
        return 0.0

    cai = math.exp(log_sum / count)
    return round(cai, 4)


def find_rare_codons(dna: str) -> list[dict[str, str | float]]:
    """Encontra codons raros na sequencia (frequencia < 10%).

    Codons raros em E. coli podem causar:
    - Pausa ribossomal
    - Frameshift
    - Truncamento prematuro
    - Formacao de corpos de inclusao

    Args:
        dna: Sequencia de DNA codificante.

    Returns:
        Lista de dicts com posicao, codon, aminoacido e frequencia.
    """
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    rare: list[dict[str, str | float]] = []

    # Construir lookup de frequencia por codon
    freq_lookup: dict[str, float] = {}
    for aa, codon_list in ECOLI_K12_CODON_TABLE.items():
        for codon, freq in codon_list:
            freq_lookup[codon] = freq

    for i, codon in enumerate(codons):
        if len(codon) < 3:
            continue
        freq = freq_lookup.get(codon, 0.0)
        if freq < RARE_CODON_THRESHOLD:
            aa = CODON_TO_AA.get(codon, "?")
            rare.append({
                "position": i + 1,
                "codon": codon,
                "amino_acid": aa,
                "frequency": freq,
            })

    return rare


def find_restriction_sites(dna: str) -> list[dict[str, str | int]]:
    """Busca sitios de restricao indesejados na sequencia completa do inserto.

    Verifica tanto os sitios de clonagem (NdeI, XhoI) internos quanto
    sitios de enzimas comuns (BamHI, EcoRI, HindIII).

    Args:
        dna: Sequencia completa do inserto (com NdeI/XhoI flanqueadores).

    Returns:
        Lista de dicts com nome da enzima, sitio e posicao.
    """
    results: list[dict[str, str | int]] = []

    # Para NdeI e XhoI, verificar apenas ocorrencias INTERNAS
    # (excluir as posicoes flanqueadoras intencionais)
    coding_region = dna[6:-6]  # Remover NdeI (6 nt) e XhoI (6 nt)

    for name, site in RESTRICTION_SITES_TO_AVOID.items():
        if name in ("NdeI", "XhoI"):
            # Buscar apenas na regiao codificante interna
            pos = coding_region.find(site)
            while pos != -1:
                results.append({
                    "enzyme": name,
                    "site": site,
                    "position": pos + 6,  # Offset pelo NdeI flanqueador
                })
                pos = coding_region.find(site, pos + 1)
        else:
            # Para outras enzimas, buscar na sequencia completa
            pos = dna.find(site)
            while pos != -1:
                results.append({
                    "enzyme": name,
                    "site": site,
                    "position": pos,
                })
                pos = dna.find(site, pos + 1)

    return results


def find_homopolymer_runs(dna: str) -> list[dict[str, str | int]]:
    """Detecta corridas homopolimericas longas (>8 nt do mesmo nucleotideo).

    Homopolimeros longos causam:
    - Erros de slippage na DNA polimerase durante sintese
    - Instabilidade da sequencia em E. coli
    - Problemas de sequenciamento

    Args:
        dna: Sequencia de DNA.

    Returns:
        Lista de dicts com posicao, nucleotideo e comprimento.
    """
    runs: list[dict[str, str | int]] = []
    pattern = re.compile(r"(.)\1{" + str(MAX_HOMOPOLYMER_LENGTH) + r",}")

    for match in pattern.finditer(dna):
        runs.append({
            "position": match.start(),
            "nucleotide": match.group(1),
            "length": len(match.group()),
        })

    return runs


def find_shine_dalgarno(dna: str) -> list[int]:
    """Busca sequencias tipo Shine-Dalgarno internas.

    AAGGAG ou variantes podem causar iniciacao espuria de traducao em E. coli,
    produzindo fragmentos proteicos truncados ou com frame incorreto.

    Args:
        dna: Sequencia de DNA codificante.

    Returns:
        Lista de posicoes onde a sequencia SD foi encontrada.
    """
    positions: list[int] = []
    pos = dna.find(SHINE_DALGARNO_PATTERN)
    while pos != -1:
        positions.append(pos)
        pos = dna.find(SHINE_DALGARNO_PATTERN, pos + 1)
    return positions


def check_rho_independent_terminator(dna: str) -> list[int]:
    """Busca potenciais terminadores rho-independentes.

    Um terminador rho-independente tipico tem:
    - Um stem-loop (palindromo GC-rico ~7-20 nt)
    - Seguido por uma corrida de Us (4+ timinas no DNA template)

    Usamos uma heuristica simplificada: buscar >= 4 Ts precedidos por
    uma regiao GC-rica de pelo menos 10 nt com >60% GC.

    Args:
        dna: Sequencia de DNA codificante.

    Returns:
        Lista de posicoes suspeitas.
    """
    suspects: list[int] = []
    # Buscar corridas de 4+ Ts (que transcrevem para poli-U)
    poly_t = re.compile(r"T{4,}")
    for match in poly_t.finditer(dna):
        t_start = match.start()
        if t_start < 10:
            continue
        # Verificar se os 20 nt anteriores sao GC-ricos
        upstream = dna[max(0, t_start - 20):t_start]
        gc_frac = sum(1 for nt in upstream if nt in "GC") / max(len(upstream), 1)
        if gc_frac > 0.60:
            suspects.append(t_start)

    return suspects


# ---------------------------------------------------------------------------
# Funcao principal de otimizacao
# ---------------------------------------------------------------------------


def run_codon_optimization(protein_sequence: str) -> CodonOptimizationResult:
    """Executa o pipeline completo de otimizacao de codons para E. coli.

    1. Traduz aminoacidos para codons de alta frequencia
    2. Remove sitios de restricao internos
    3. Adiciona elementos flanqueadores (NdeI, stop, XhoI)
    4. Calcula metricas de qualidade (CAI, GC%, codons raros)
    5. Verifica elementos problematicos

    Args:
        protein_sequence: Sequencia proteica completa do constructo.

    Returns:
        CodonOptimizationResult com DNA otimizado e metricas.
    """
    # Passo 1-2: Otimizar codons e remover sitios de restricao
    coding_dna = optimize_codons(protein_sequence)

    # Passo 3: Adicionar elementos flanqueadores
    full_insert = add_flanking_elements(coding_dna)

    # Passo 4: Metricas de qualidade sobre a regiao codificante
    gc = calculate_gc_content(coding_dna)
    cai = calculate_cai(coding_dna)
    rare = find_rare_codons(coding_dna)

    # Passo 5: Verificacoes sobre o inserto completo
    sites = find_restriction_sites(full_insert)
    homopolymers = find_homopolymer_runs(full_insert)
    sd_positions = find_shine_dalgarno(coding_dna)
    rho_positions = check_rho_independent_terminator(coding_dna)

    # Classificar sitios NdeI/XhoI internos
    has_ndei = any(s["enzyme"] == "NdeI" for s in sites)
    has_xhoi = any(s["enzyme"] == "XhoI" for s in sites)

    # Gerar avisos
    warnings: list[str] = []
    if gc < 0.40 or gc > 0.60:
        warnings.append(
            f"GC content = {gc:.1%} fora do range ideal (40-60%)"
        )
    if cai < 0.80:
        warnings.append(
            f"CAI = {cai:.3f} abaixo do ideal (>0.80)"
        )
    if rare:
        warnings.append(
            f"{len(rare)} codon(s) raro(s) detectado(s)"
        )
    if has_ndei:
        warnings.append("Sitio NdeI interno detectado na regiao codificante")
    if has_xhoi:
        warnings.append("Sitio XhoI interno detectado na regiao codificante")
    if sd_positions:
        warnings.append(
            f"{len(sd_positions)} sequencia(s) Shine-Dalgarno interna(s) encontrada(s)"
        )
    if rho_positions:
        warnings.append(
            f"{len(rho_positions)} potencial(is) terminador(es) rho-independente(s)"
        )
    if homopolymers:
        warnings.append(
            f"{len(homopolymers)} corrida(s) homopolimerica(s) >8 nt"
        )

    return CodonOptimizationResult(
        dna_sequence=full_insert,
        length_nt=len(full_insert),
        gc_content=gc,
        cai=cai,
        rare_codons=rare,
        restriction_sites_found=sites,
        has_internal_ndei=has_ndei,
        has_internal_xhoi=has_xhoi,
        has_shine_dalgarno=len(sd_positions) > 0,
        homopolymer_runs=homopolymers,
        warnings=warnings,
    )


def to_dict(result: CodonOptimizationResult) -> dict:
    """Serializa o resultado da otimizacao para JSON.

    Args:
        result: Resultado da otimizacao de codons.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "dna_sequence": result.dna_sequence,
        "length_nt": result.length_nt,
        "gc_content": result.gc_content,
        "cai": result.cai,
        "rare_codons": result.rare_codons,
        "restriction_sites_found": result.restriction_sites_found,
        "has_internal_ndei": result.has_internal_ndei,
        "has_internal_xhoi": result.has_internal_xhoi,
        "has_shine_dalgarno": result.has_shine_dalgarno,
        "homopolymer_runs": result.homopolymer_runs,
        "warnings": result.warnings,
    }
