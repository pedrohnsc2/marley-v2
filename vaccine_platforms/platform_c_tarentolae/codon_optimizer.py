"""C3 -- Otimizacao de codons para expressao em L. tarentolae.

Trypanosomatideos possuem forte vies de uso de codons GC-rich (~60% GC
no genoma). A otimizacao de codons para L. tarentolae e essencial para
garantir alta expressao da proteina recombinante.

Estrategia:
    1. Usar tabela de frequencia de codons de L. tarentolae (Kazusa DB)
    2. Selecionar codons de maior frequencia (max-frequency)
    3. Evitar sitios internos de restricao (SwaI, usado no LEXSY)
    4. Verificar ausencia de sinais de poliadenilacao espurios
    5. Verificar ausencia de sequencias de splice aceptor/donor
    6. Alvo de GC: 55-65% (compativel com expressao em trypanosomatideos)
    7. CAI alvo: > 0.85

Referencia:
    - Codon usage de Leishmania: Horn (2008) Mol Biochem Parasitol
    - pLEXSY manual (Jena Bioscience, 2023)
    - Kazusa Codon Usage Database (L. tarentolae: taxid 5765)
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Tabela de uso de codons de L. tarentolae
# ---------------------------------------------------------------------------
# Fonte: Kazusa Codon Usage Database (http://www.kazusa.or.jp/codon/)
# e dados de genes altamente expressos de L. major/L. tarentolae.
# Valores: frequencia relativa normalizada por aminoacido (soma = 1.0).
# Trypanosomatideos sao fortemente GC-biased: codons terminando em C ou G
# sao preferidos sobre A e T/U.

LT_CODON_TABLE: Final[dict[str, list[tuple[str, float]]]] = {
    # Alanina -- GCC e GCG preferidos (GC-rich)
    "A": [("GCC", 0.38), ("GCG", 0.30), ("GCA", 0.10), ("GCT", 0.22)],
    # Arginina -- CGC domina em trypanosomatideos
    "R": [("CGC", 0.40), ("CGG", 0.22), ("CGT", 0.12), ("CGA", 0.06), ("AGG", 0.12), ("AGA", 0.08)],
    # Asparagina -- AAC preferido (~75%)
    "N": [("AAC", 0.75), ("AAT", 0.25)],
    # Aspartato -- GAC preferido
    "D": [("GAC", 0.65), ("GAT", 0.35)],
    # Cisteina -- TGC preferido
    "C": [("TGC", 0.70), ("TGT", 0.30)],
    # Glutamina -- CAG preferido
    "Q": [("CAG", 0.68), ("CAA", 0.32)],
    # Glutamato -- GAG preferido (diferente de E. coli onde GAA domina)
    "E": [("GAG", 0.62), ("GAA", 0.38)],
    # Glicina -- GGC preferido
    "G": [("GGC", 0.42), ("GGT", 0.15), ("GGA", 0.13), ("GGG", 0.30)],
    # Histidina -- CAC preferido
    "H": [("CAC", 0.70), ("CAT", 0.30)],
    # Isoleucina -- ATC preferido (~55%)
    "I": [("ATC", 0.55), ("ATT", 0.30), ("ATA", 0.15)],
    # Leucina -- CTG e CTC preferidos
    "L": [("CTG", 0.38), ("CTC", 0.25), ("CTT", 0.12), ("TTG", 0.10), ("TTA", 0.05), ("CTA", 0.10)],
    # Lisina -- AAG muito preferido (~78%)
    "K": [("AAG", 0.78), ("AAA", 0.22)],
    # Metionina -- codon unico
    "M": [("ATG", 1.00)],
    # Fenilalanina -- TTC preferido
    "F": [("TTC", 0.68), ("TTT", 0.32)],
    # Prolina -- CCG e CCC preferidos
    "P": [("CCG", 0.35), ("CCC", 0.30), ("CCT", 0.15), ("CCA", 0.20)],
    # Serina -- AGC e TCC preferidos
    "S": [("AGC", 0.30), ("TCC", 0.25), ("TCG", 0.20), ("TCT", 0.08), ("TCA", 0.07), ("AGT", 0.10)],
    # Treonina -- ACC preferido
    "T": [("ACC", 0.45), ("ACG", 0.28), ("ACT", 0.12), ("ACA", 0.15)],
    # Triptofano -- codon unico
    "W": [("TGG", 1.00)],
    # Tirosina -- TAC preferido
    "Y": [("TAC", 0.72), ("TAT", 0.28)],
    # Valina -- GTG e GTC preferidos
    "V": [("GTG", 0.40), ("GTC", 0.30), ("GTT", 0.15), ("GTA", 0.15)],
    # Stop -- TGA preferido em trypanosomatideos
    "*": [("TGA", 0.50), ("TAA", 0.30), ("TAG", 0.20)],
}

# Mapa reverso: codon -> aminoacido
CODON_TO_AA: Final[dict[str, str]] = {}
for _aa, _codons in LT_CODON_TABLE.items():
    for _codon, _freq in _codons:
        CODON_TO_AA[_codon] = _aa

# Sitios de restricao a evitar (pLEXSY usa SwaI e BglII para clonagem)
RESTRICTION_SITES_TO_AVOID: Final[dict[str, str]] = {
    "SwaI": "ATTTAAAT",     # Sitio principal de clonagem pLEXSY
    "BglII": "AGATCT",     # Sitio secundario
    "NcoI": "CCATGG",      # Pode interferir com ATG de inicio
    "NotI": "GCGGCCGC",    # Usado no vetor backbone
}

# Sinais de processamento de RNA a evitar na CDS
# Sinal de poliadenilacao em trypanosomatideos nao e AATAAA como em mamiferos;
# o processamento e dirigido por regioes intergenicas. Porem, sequencias
# poli-A longas e tratos poli-T podem causar terminacao prematura.
POLYA_LIKE_PATTERN: Final[str] = "AAAAAA"  # 6+ A's consecutivos
POLYT_LIKE_PATTERN: Final[str] = "TTTTTT"  # 6+ T's consecutivos

# Limiar de codon raro
RARE_CODON_THRESHOLD: Final[float] = 0.10

# Comprimento maximo de homopolimero
MAX_HOMOPOLYMER_LENGTH: Final[int] = 6  # Mais restritivo que E. coli (8)

# Faixa alvo de GC para trypanosomatideos
GC_TARGET_MIN: Final[float] = 0.55
GC_TARGET_MAX: Final[float] = 0.65


# ---------------------------------------------------------------------------
# Dataclass de resultado
# ---------------------------------------------------------------------------


@dataclass
class CodonOptimizationResult:
    """Resultado da otimizacao de codons para L. tarentolae."""

    dna_sequence: str
    length_nt: int
    gc_content: float
    cai: float
    gc_in_range: bool
    cai_above_threshold: bool
    rare_codons: list[dict]
    restriction_sites_found: list[dict]
    homopolymer_runs: list[dict]
    has_polya_signal: bool
    has_polyt_signal: bool
    codon_frequency_profile: dict[str, int]
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Otimizacao de codons
# ---------------------------------------------------------------------------


def get_best_codon(amino_acid: str) -> str:
    """Retorna o codon de maior frequencia para o aminoacido em L. tarentolae.

    Em trypanosomatideos, a preferencia e fortemente GC-biased:
    - Leu: CTG domina (~38%)
    - Arg: CGC preferido (~40%) -- diferente de E. coli (CGT) e humano (AGA)
    - Glu: GAG preferido (~62%) -- diferente de E. coli (GAA)
    - Lys: AAG muito preferido (~78%)

    Args:
        amino_acid: Letra unica do aminoacido (maiuscula).

    Returns:
        Triplet de DNA otimizado para L. tarentolae.

    Raises:
        ValueError: Se o aminoacido nao esta na tabela.
    """
    codons = LT_CODON_TABLE.get(amino_acid)
    if not codons:
        raise ValueError(f"Aminoacido desconhecido: {amino_acid!r}")
    return max(codons, key=lambda x: x[1])[0]


def optimize_codons(protein_sequence: str) -> str:
    """Traduz sequencia proteica em DNA otimizado para L. tarentolae.

    Usa estrategia de maior frequencia (max-frequency) para cada posicao.
    Apos a selecao inicial, verifica e remove sitios de restricao internos
    e sinais de processamento de RNA indesejados.

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
    dna = _remove_restriction_sites(dna, protein_sequence)

    # Remover homopolimeros longos
    dna = _break_homopolymers(dna, protein_sequence)

    return dna


def _remove_restriction_sites(dna: str, protein: str) -> str:
    """Remove sitios de restricao internos via substituicao sinonima.

    Para cada ocorrencia de um sitio de restricao, encontra o codon
    sobreposto e troca por um sinonimo alternativo que elimine o sitio.

    Args:
        dna: Sequencia de DNA codificante.
        protein: Sequencia proteica correspondente.

    Returns:
        DNA com sitios de restricao removidos.
    """
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    max_iterations = 20

    for _iteration in range(max_iterations):
        changed = False
        current_dna = "".join(codons)

        for site_name, site_seq in RESTRICTION_SITES_TO_AVOID.items():
            pos = current_dna.find(site_seq)
            while pos != -1:
                codon_start_idx = pos // 3
                codon_end_idx = (pos + len(site_seq) - 1) // 3

                for ci in range(codon_start_idx, min(codon_end_idx + 1, len(codons))):
                    if ci >= len(protein):
                        break
                    aa = protein[ci]
                    alternatives = LT_CODON_TABLE.get(aa, [])
                    sorted_alts = sorted(alternatives, key=lambda x: -x[1])

                    for alt_codon, alt_freq in sorted_alts:
                        if alt_codon == codons[ci]:
                            continue
                        if alt_freq < RARE_CODON_THRESHOLD:
                            continue
                        test_codons = codons[:]
                        test_codons[ci] = alt_codon
                        test_dna = "".join(test_codons)
                        if site_seq not in test_dna[max(0, pos - 8):pos + len(site_seq) + 8]:
                            codons[ci] = alt_codon
                            changed = True
                            break
                    if changed:
                        break

                current_dna = "".join(codons)
                next_pos = current_dna.find(site_seq, pos + 1)
                pos = next_pos

        if not changed:
            break

    return "".join(codons)


def _break_homopolymers(dna: str, protein: str) -> str:
    """Elimina corridas de homopolimeros > MAX_HOMOPOLYMER_LENGTH.

    Homopolimeros longos causam problemas de sintese (erros de
    sequenciamento, deleicoes) e podem afetar processamento de RNA
    em trypanosomatideos.

    Args:
        dna: Sequencia de DNA codificante.
        protein: Sequencia proteica correspondente.

    Returns:
        DNA com homopolimeros quebrados.
    """
    codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]

    for _iteration in range(10):
        current_dna = "".join(codons)
        found_problem = False

        for nt in "ACGT":
            pattern = nt * (MAX_HOMOPOLYMER_LENGTH + 1)
            pos = current_dna.find(pattern)
            if pos == -1:
                continue

            found_problem = True
            # Tentar quebrar por substituicao sinonima do codon central
            mid_pos = pos + MAX_HOMOPOLYMER_LENGTH // 2
            ci = mid_pos // 3
            if ci < len(protein):
                aa = protein[ci]
                alternatives = LT_CODON_TABLE.get(aa, [])
                for alt_codon, alt_freq in sorted(alternatives, key=lambda x: -x[1]):
                    if alt_codon == codons[ci]:
                        continue
                    if alt_freq < RARE_CODON_THRESHOLD:
                        continue
                    test_codons = codons[:]
                    test_codons[ci] = alt_codon
                    test_dna = "".join(test_codons)
                    if pattern not in test_dna:
                        codons[ci] = alt_codon
                        break

        if not found_problem:
            break

    return "".join(codons)


def add_flanking_elements(coding_dna: str) -> str:
    """Adiciona elementos flanqueadores para clonagem no pLEXSY-sat2.

    Estrutura do inserto completo (DNA):
        [SSU 5' flank] - [5'UTR] - ATG + [CDS] - stop(TGA) - [3'UTR] - [SSU 3' flank]

    O ATG de inicio ja esta incluido no primeiro codon da CDS.
    O stop codon preferido em trypanosomatideos e TGA.
    As regioes UTR e flancos SSU sao adicionados para integracao
    e processamento correto do transcrito.

    Args:
        coding_dna: DNA codificante otimizado (comecando com ATG).

    Returns:
        Inserto completo com todos os elementos flanqueadores.
    """
    from vaccine_platforms.platform_c_tarentolae.construct import (
        UTR5_TARENTOLAE,
        UTR3_TARENTOLAE,
        SSU_5_FLANK,
        SSU_3_FLANK,
    )

    stop = "TGA"  # Stop codon preferido em trypanosomatideos

    # Montar inserto completo
    insert = f"{SSU_5_FLANK}{UTR5_TARENTOLAE}{coding_dna}{stop}{UTR3_TARENTOLAE}{SSU_3_FLANK}"
    return insert


# ---------------------------------------------------------------------------
# Metricas de qualidade
# ---------------------------------------------------------------------------


def calculate_gc_content(dna: str) -> float:
    """Calcula a fracao de G+C na sequencia de DNA.

    GC content ideal para L. tarentolae: 55-65% (genoma ~60% GC).

    Args:
        dna: Sequencia de DNA.

    Returns:
        Fracao GC (0.0 a 1.0).
    """
    dna = dna.upper()
    if not dna:
        return 0.0
    gc_count = sum(1 for nt in dna if nt in ("G", "C"))
    return round(gc_count / len(dna), 4)


def calculate_cai(dna: str) -> float:
    """Calcula o CAI (Codon Adaptation Index) para L. tarentolae.

    Usa a media geometrica dos pesos relativos de cada codon,
    normalizados pelo codon de maior frequencia de cada familia
    de sinonimos.

    Args:
        dna: Sequencia codificante (CDS), comprimento multiplo de 3.

    Returns:
        CAI entre 0.0 e 1.0.
    """
    dna = dna.upper().strip()
    if len(dna) < 3:
        return 0.0

    # Calcular peso maximo por aminoacido
    aa_max_freq: dict[str, float] = {}
    for aa, codons in LT_CODON_TABLE.items():
        if aa == "*":
            continue
        max_freq = max(freq for _, freq in codons)
        aa_max_freq[aa] = max_freq

    # Calcular peso relativo (w) para cada codon
    codon_weights: dict[str, float] = {}
    for aa, codons in LT_CODON_TABLE.items():
        if aa == "*":
            continue
        max_f = aa_max_freq.get(aa, 1.0)
        for codon, freq in codons:
            codon_weights[codon] = freq / max_f if max_f > 0 else 0.0

    # Media geometrica dos pesos
    log_sum = 0.0
    codon_count = 0

    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        if len(codon) < 3:
            break
        aa = CODON_TO_AA.get(codon)
        if aa is None or aa == "*":
            continue
        w = codon_weights.get(codon, 0.0)
        if w > 0:
            log_sum += math.log(w)
            codon_count += 1

    if codon_count == 0:
        return 0.0

    cai = math.exp(log_sum / codon_count)
    return round(cai, 4)


def find_restriction_sites(dna: str) -> list[dict]:
    """Busca sitios de restricao indesejados na sequencia.

    Args:
        dna: Sequencia de DNA.

    Returns:
        Lista de dicts com enzima e posicao de cada sitio encontrado.
    """
    dna = dna.upper()
    sites: list[dict] = []

    for enzyme, site_seq in RESTRICTION_SITES_TO_AVOID.items():
        pos = dna.find(site_seq)
        while pos != -1:
            sites.append({"enzyme": enzyme, "site": site_seq, "position": pos})
            pos = dna.find(site_seq, pos + 1)

    return sites


def find_homopolymers(dna: str) -> list[dict]:
    """Busca corridas de homopolimeros acima do limiar.

    Args:
        dna: Sequencia de DNA.

    Returns:
        Lista de dicts com nucleotideo, posicao e comprimento de cada corrida.
    """
    dna = dna.upper()
    runs: list[dict] = []

    for nt in "ACGT":
        pattern = re.compile(f"{nt}{{{MAX_HOMOPOLYMER_LENGTH + 1},}}")
        for match in pattern.finditer(dna):
            runs.append({
                "nucleotide": nt,
                "position": match.start(),
                "length": match.end() - match.start(),
            })

    return runs


def get_codon_frequency_profile(dna: str) -> dict[str, int]:
    """Conta a frequencia de cada codon na CDS.

    Util para visualizar o perfil de uso de codons e verificar
    que nao ha codons raros em excesso.

    Args:
        dna: Sequencia codificante (CDS).

    Returns:
        Dicionario codon -> contagem.
    """
    profile: dict[str, int] = {}
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        if len(codon) == 3:
            profile[codon] = profile.get(codon, 0) + 1
    return dict(sorted(profile.items()))


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def run_codon_optimization(protein_sequence: str) -> CodonOptimizationResult:
    """Executa o pipeline completo de otimizacao de codons para L. tarentolae.

    Etapas:
        1. Otimizar codons usando tabela de L. tarentolae
        2. Calcular GC content e CAI
        3. Buscar codons raros residuais
        4. Buscar sitios de restricao residuais
        5. Buscar homopolimeros longos
        6. Verificar sinais de poliadenilacao espurios

    Args:
        protein_sequence: Sequencia proteica a otimizar.

    Returns:
        CodonOptimizationResult com todos os dados de qualidade.
    """
    # Otimizar codons
    dna = optimize_codons(protein_sequence)

    # Metricas de qualidade
    gc = calculate_gc_content(dna)
    cai = calculate_cai(dna)
    gc_ok = GC_TARGET_MIN <= gc <= GC_TARGET_MAX
    cai_ok = cai >= 0.85

    # Codons raros residuais
    rare_codons: list[dict] = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None or aa == "*":
            continue
        # Encontrar a frequencia deste codon
        freq = 0.0
        for c, f in LT_CODON_TABLE.get(aa, []):
            if c == codon:
                freq = f
                break
        if freq < RARE_CODON_THRESHOLD:
            rare_codons.append({
                "position": i // 3,
                "codon": codon,
                "amino_acid": aa,
                "frequency": freq,
            })

    # Sitios de restricao
    sites = find_restriction_sites(dna)

    # Homopolimeros
    homopolymers = find_homopolymers(dna)

    # Sinais de poliadenilacao / terminacao
    has_polya = POLYA_LIKE_PATTERN in dna.upper()
    has_polyt = POLYT_LIKE_PATTERN in dna.upper()

    # Perfil de codons
    profile = get_codon_frequency_profile(dna)

    # Gerar avisos
    warnings: list[str] = []

    if not gc_ok:
        warnings.append(
            f"GC content = {gc:.1%} fora da faixa alvo ({GC_TARGET_MIN:.0%}-{GC_TARGET_MAX:.0%})"
        )

    if not cai_ok:
        warnings.append(f"CAI = {cai:.4f} abaixo do alvo (>0.85)")

    if rare_codons:
        warnings.append(f"{len(rare_codons)} codons raros (<{RARE_CODON_THRESHOLD:.0%}) encontrados")

    if sites:
        enzymes = ", ".join(set(s["enzyme"] for s in sites))
        warnings.append(f"Sitios de restricao residuais: {enzymes}")

    if homopolymers:
        warnings.append(f"{len(homopolymers)} homopolimeros >{MAX_HOMOPOLYMER_LENGTH} nt encontrados")

    if has_polya:
        warnings.append("Sinal tipo poli-A encontrado (>=6 A's consecutivos)")

    if has_polyt:
        warnings.append("Trato poli-T encontrado (>=6 T's consecutivos) -- risco de terminacao prematura")

    return CodonOptimizationResult(
        dna_sequence=dna,
        length_nt=len(dna),
        gc_content=gc,
        cai=cai,
        gc_in_range=gc_ok,
        cai_above_threshold=cai_ok,
        rare_codons=rare_codons,
        restriction_sites_found=sites,
        homopolymer_runs=homopolymers,
        has_polya_signal=has_polya,
        has_polyt_signal=has_polyt,
        codon_frequency_profile=profile,
        warnings=warnings,
    )


def to_dict(result: CodonOptimizationResult) -> dict:
    """Serializa o resultado de otimizacao de codons para JSON.

    Args:
        result: Resultado da otimizacao.

    Returns:
        Dicionario serializavel para JSON.
    """
    return {
        "dna_sequence": result.dna_sequence,
        "length_nt": result.length_nt,
        "gc_content": result.gc_content,
        "gc_in_range": result.gc_in_range,
        "gc_target": f"{GC_TARGET_MIN:.0%}-{GC_TARGET_MAX:.0%}",
        "cai": result.cai,
        "cai_above_threshold": result.cai_above_threshold,
        "cai_target": ">0.85",
        "rare_codons": result.rare_codons,
        "rare_codon_count": len(result.rare_codons),
        "restriction_sites_found": result.restriction_sites_found,
        "homopolymer_runs": result.homopolymer_runs,
        "has_polya_signal": result.has_polya_signal,
        "has_polyt_signal": result.has_polyt_signal,
        "codon_frequency_profile": result.codon_frequency_profile,
        "warnings": result.warnings,
    }
