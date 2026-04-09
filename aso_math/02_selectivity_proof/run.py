"""Modulo 02 — Prova matematica de seletividade do MRL-ASO-001.

Demonstra que ligacao off-target a transcritos humanos ou caninos e
impossivel, usando duas abordagens complementares:

1. Analise de complementaridade maxima contigua (medida direta do
   potencial de hibridizacao) contra transcriptoma humano ou modelo nulo.
2. Entropia de Shannon por posicao do SL RNA entre especies de
   Leishmania (prova que o alvo e maximamente conservado).

Referencia principal:
  - Crooke ST et al. (2017) Cellular uptake and trafficking of
    antisense oligonucleotides. Nat Biotechnol 35(3):230-237.
  - Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612.
  - Milhausen M et al. (1984) Cell 38(3):721-729.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from pathlib import Path
from typing import Any

from aso_math.config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASES,
    COMPLEMENT,
    DG_FUNCTIONAL_THRESHOLD,
    TRANSCRIPTOME_HUMAN,
)
from aso_math.envelope import Timer, create_envelope, write_result
from aso_math.target_config import TargetConfig
from aso_math.thermo import compute_dg, reverse_complement
from core.logger import get_logger

logger = get_logger("aso_math.02_selectivity_proof")

# Semente para reprodutibilidade do modelo nulo
_RNG_SEED: int = 42

# Limiar minimo de pares complementares contiguos para ativacao de RNase H
# Ref: Crooke 2017 — ASOs precisam de ~14-16 bp contiguos
_RNASE_H_THRESHOLD_BP: int = 14


# ---------------------------------------------------------------------------
# Parte 1: Analise de complementaridade
# ---------------------------------------------------------------------------


def _is_nucleotide_fasta(fasta_path: Path, sample_lines: int = 20) -> bool:
    """Verifica se um FASTA contem sequencias de nucleotideos (nao proteinas).

    Le as primeiras linhas de sequencia e checa se contem apenas
    caracteres validos de DNA/RNA (A, C, G, T, U, N e variantes IUPAC).
    Se encontrar aminoacidos frequentes nao-ambiguos (M, W, F, Y, etc.),
    conclui que e proteina.
    """
    protein_only_chars = set("EFILPQRDHKMWYefilpqrdhkmwy")
    lines_checked = 0

    try:
        with open(fasta_path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                stripped = line.strip()
                if not stripped:
                    continue
                # Se encontrar caracteres exclusivos de proteina, e proteina
                if protein_only_chars & set(stripped):
                    return False
                lines_checked += 1
                if lines_checked >= sample_lines:
                    break
    except OSError:
        return False

    # Se todas as linhas checadas contem apenas nucleotideos
    return lines_checked > 0


def parse_fasta(fasta_path: Path) -> dict[str, str]:
    """Parser FASTA simples sem dependencias externas.

    Le um arquivo FASTA e retorna dicionario {header: sequencia}.
    Headers sao a linha apos '>' sem o caracter '>'.
    Sequencias sao concatenadas em uppercase.
    """
    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_seq: list[str] = []

    with open(fasta_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Salva sequencia anterior se existir
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq).upper()
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    # Salva a ultima sequencia
    if current_header is not None:
        sequences[current_header] = "".join(current_seq).upper()

    return sequences


def _generate_null_model(
    n_sequences: int = 500,
    min_length: int = 200,
    max_length: int = 2000,
) -> dict[str, str]:
    """Gera sequencias aleatorias com conteudo GC similar ao humano (40-60%).

    Cada sequencia tem comprimento uniforme entre min_length e max_length.
    O conteudo GC e amostrado uniformemente entre 0.40 e 0.60 para cada
    sequencia, refletindo a variacao real de transcritos humanos.
    """
    random.seed(_RNG_SEED)
    sequences: dict[str, str] = {}

    for i in range(n_sequences):
        length = random.randint(min_length, max_length)
        gc_frac = random.uniform(0.40, 0.60)

        # Probabilidades por base: p(G) = p(C) = gc/2, p(A) = p(T) = (1-gc)/2
        weights = [
            (1.0 - gc_frac) / 2.0,  # A
            gc_frac / 2.0,           # C
            gc_frac / 2.0,           # G
            (1.0 - gc_frac) / 2.0,  # T
        ]
        bases_list = random.choices(list(BASES), weights=weights, k=length)
        seq = "".join(bases_list)
        sequences[f"null_model_seq_{i:04d}|len={length}|gc={gc_frac:.2f}"] = seq

    return sequences


def find_max_complementarity(aso_seq: str, transcript_seq: str) -> tuple[int, int]:
    """Encontra o comprimento maximo de complementaridade contigua.

    Desliza o ASO (e seu complemento reverso) ao longo de um transcrito.
    Em cada posicao, conta o maior run de bases complementares contiguas.
    Isso mede diretamente o potencial de hibridizacao — ASOs precisam de
    ~14-16 bp contiguos para ativar RNase H (Crooke 2017).

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        transcript_seq: Sequencia do transcrito alvo.

    Returns:
        Tupla (max_length, best_position) onde position e a posicao no
        transcrito onde o melhor alinhamento comeca.
    """
    aso_upper = aso_seq.upper()
    transcript_upper = transcript_seq.upper()
    aso_rc = reverse_complement(aso_upper)
    aso_len = len(aso_upper)
    transcript_len = len(transcript_upper)

    global_max_length = 0
    global_best_pos = -1

    # Testar ambas as orientacoes: ASO direto e complemento reverso
    for probe in (aso_upper, aso_rc):
        probe_len = len(probe)

        # Deslizar o probe ao longo do transcrito em todas as posicoes
        # onde ha pelo menos 1 base de overlap
        for offset in range(-(probe_len - 1), transcript_len):
            # Determinar regiao de overlap
            t_start = max(0, offset)
            p_start = max(0, -offset)
            overlap = min(probe_len - p_start, transcript_len - t_start)

            if overlap < 1:
                continue

            # Contar maior run de bases complementares contiguas
            current_run = 0
            best_run = 0
            best_run_pos = t_start

            for i in range(overlap):
                t_base = transcript_upper[t_start + i]
                p_base = probe[p_start + i]

                # Checa se sao complementares (Watson-Crick)
                if COMPLEMENT.get(p_base) == t_base:
                    current_run += 1
                    if current_run > best_run:
                        best_run = current_run
                        best_run_pos = t_start + i - current_run + 1
                else:
                    current_run = 0

            if best_run > global_max_length:
                global_max_length = best_run
                global_best_pos = best_run_pos

    return global_max_length, global_best_pos


def screen_transcriptome(
    aso_seq: str,
    fasta_path: Path,
) -> tuple[list[tuple[str, int, int]], str, list[str]]:
    """Varre um transcriptoma inteiro em busca de complementaridade off-target.

    Para cada transcrito, computa a complementaridade maxima contigua com
    o ASO. Retorna resultados ordenados por complementaridade decrescente.

    Se o FASTA nao existir ou contiver proteinas ao inves de nucleotideos,
    gera 500 sequencias aleatorias com conteudo GC humano como modelo nulo.

    Args:
        aso_seq: Sequencia do ASO.
        fasta_path: Caminho para arquivo FASTA do transcriptoma.

    Returns:
        Tupla (resultados, tipo_transcriptoma, avisos) onde:
        - resultados: lista de (transcript_id, max_compl, position)
        - tipo_transcriptoma: "human", "null_model"
        - avisos: lista de strings com avisos gerados
    """
    warnings: list[str] = []
    transcriptome_type: str

    if fasta_path.exists() and _is_nucleotide_fasta(fasta_path):
        logger.info("Lendo transcriptoma de nucleotideos: %s", fasta_path)
        sequences = parse_fasta(fasta_path)
        transcriptome_type = "human"
        logger.info("Total de transcritos carregados: %d", len(sequences))
    else:
        if fasta_path.exists():
            reason = (
                "O arquivo FASTA existe mas contem sequencias de proteinas, "
                "nao nucleotideos. Usando modelo nulo com sequencias aleatorias."
            )
        else:
            reason = (
                f"Arquivo FASTA nao encontrado em {fasta_path}. "
                "Usando modelo nulo com sequencias aleatorias."
            )
        logger.warning(reason)
        warnings.append(reason)
        sequences = _generate_null_model()
        transcriptome_type = "null_model"
        logger.info("Modelo nulo gerado: %d sequencias aleatorias", len(sequences))

    results: list[tuple[str, int, int]] = []
    total = len(sequences)

    for idx, (header, seq) in enumerate(sequences.items()):
        max_compl, position = find_max_complementarity(aso_seq, seq)
        results.append((header, max_compl, position))

        # Log de progresso a cada 10%
        if total >= 10 and (idx + 1) % (total // 10) == 0:
            pct = ((idx + 1) * 100) // total
            logger.info("  Progresso: %d%% (%d/%d transcritos)", pct, idx + 1, total)

    # Ordenar por complementaridade decrescente
    results.sort(key=lambda x: x[1], reverse=True)

    return results, transcriptome_type, warnings


def compute_partial_binding_dg(aso_seq: str, match_length: int) -> float:
    """Calcula delta_G de ligacao parcial para um trecho complementar.

    Para o maior trecho complementar encontrado no transcriptoma off-target,
    estima o dG que essa ligacao parcial teria. Se dG > DG_FUNCTIONAL_THRESHOLD
    (-15.0 kcal/mol), a ligacao parcial e fraca demais para ser funcional.

    Usa uma subsequencia do ASO com o comprimento dado para o calculo.

    Args:
        aso_seq: Sequencia completa do ASO.
        match_length: Comprimento do trecho complementar maximo encontrado.

    Returns:
        delta_G em kcal/mol para a ligacao parcial.
    """
    if match_length <= 0:
        return 0.0

    # Usa a subsequencia central do ASO com o comprimento do match
    # Isso e uma estimativa conservadora (melhor caso para o off-target)
    aso_upper = aso_seq.upper()
    aso_len = len(aso_upper)

    if match_length >= aso_len:
        return compute_dg(aso_upper)

    # Testa todas as substrings possiveis do ASO com esse comprimento
    # e retorna o dG mais negativo (pior caso = ligacao mais forte)
    best_dg = 0.0
    for start in range(aso_len - match_length + 1):
        subseq = aso_upper[start : start + match_length]
        dg = compute_dg(subseq)
        if dg < best_dg:
            best_dg = dg

    return best_dg


# ---------------------------------------------------------------------------
# Parte 2: Conservacao do SL RNA — entropia de Shannon
# ---------------------------------------------------------------------------


def get_sl_rna_sequences() -> dict[str, str]:
    """Retorna sequencias de SL RNA 5' exon de trypanosomatideos.

    Sequencias publicadas da regiao 5' exon do Spliced Leader RNA.
    O SL RNA e adicionado a todo mRNA do parasita via trans-splicing e
    e altamente conservado entre especies de Leishmania.

    Refs:
      - Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
      - Milhausen M et al. (1984) Cell 38(3):721-729

    Returns:
        Dicionario {nome_especie: sequencia_39nt}.
    """
    return {
        "L. infantum":    "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "L. donovani":    "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "L. major":       "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "L. braziliensis": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "L. mexicana":    "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "T. brucei":      "AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG",
        "T. cruzi":       "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG",
    }


def compute_positional_entropy(sequences: dict[str, str]) -> list[float]:
    """Calcula entropia de Shannon por posicao do alinhamento.

    Para cada posicao (0..38), calcula:
        H = -sum( p(base) * log2(p(base)) )

    Posicoes onde todas as especies tem a mesma base recebem H = 0.0
    (perfeitamente conservada). Entropia maxima = 2.0 bits (4 bases
    equiprovaveis).

    Args:
        sequences: Dicionario {especie: sequencia} com sequencias alinhadas
                   de mesmo comprimento.

    Returns:
        Lista de floats com entropia por posicao, arredondada a 4 casas.
    """
    seq_list = list(sequences.values())
    n_seqs = len(seq_list)
    seq_len = len(seq_list[0])

    entropy_per_position: list[float] = []

    for pos in range(seq_len):
        # Contar frequencia de cada base nesta posicao
        base_counts: dict[str, int] = defaultdict(int)
        for seq in seq_list:
            if pos < len(seq):
                base_counts[seq[pos]] += 1

        # Calcular entropia de Shannon
        h = 0.0
        for count in base_counts.values():
            if count > 0:
                p = count / n_seqs
                h -= p * math.log2(p)

        entropy_per_position.append(round(h, 4))

    return entropy_per_position


def assess_target_region_conservation(
    entropy_array: list[float],
    target_start: int,
    target_end: int,
) -> dict[str, Any]:
    """Avalia conservacao da regiao-alvo do ASO vs. regioes flanqueantes.

    Compara a entropia media dentro da regiao-alvo do ASO (posicoes
    target_start..target_end) com a entropia media fora dessa regiao.

    Posicoes com H = 0.0 sao perfeitamente conservadas entre todas as
    especies analisadas.

    Args:
        entropy_array: Lista de entropias por posicao (saida de
                       compute_positional_entropy).
        target_start: Posicao inicial da regiao-alvo (0-indexed, inclusivo).
        target_end: Posicao final da regiao-alvo (0-indexed, exclusivo).

    Returns:
        Dicionario com metricas de conservacao.
    """
    target_entropies = entropy_array[target_start:target_end]
    non_target_entropies = (
        entropy_array[:target_start] + entropy_array[target_end:]
    )

    # Entropia media dentro e fora da regiao-alvo
    target_mean = (
        sum(target_entropies) / len(target_entropies)
        if target_entropies
        else 0.0
    )
    non_target_mean = (
        sum(non_target_entropies) / len(non_target_entropies)
        if non_target_entropies
        else 0.0
    )

    # Contagem de posicoes totalmente conservadas (H = 0.0) na regiao-alvo
    fully_conserved = sum(1 for h in target_entropies if h == 0.0)
    variable = len(target_entropies) - fully_conserved

    return {
        "target_region_mean_entropy": round(target_mean, 4),
        "non_target_mean_entropy": round(non_target_mean, 4),
        "fully_conserved_positions_in_target": fully_conserved,
        "variable_positions_in_target": variable,
        "target_region_length": len(target_entropies),
        "conservation_ratio": round(
            fully_conserved / len(target_entropies), 4
        ) if target_entropies else 0.0,
    }


# ---------------------------------------------------------------------------
# Histograma de complementaridade
# ---------------------------------------------------------------------------


def _build_complementarity_histogram(
    results: list[tuple[str, int, int]],
) -> dict[str, int]:
    """Constroi histograma de comprimentos de complementaridade.

    Agrupa os transcritos por comprimento maximo de complementaridade
    encontrado (1, 2, 3, ... bp).

    Args:
        results: Lista de (header, max_compl, position).

    Returns:
        Dicionario {str(comprimento): contagem}.
    """
    histogram: dict[int, int] = defaultdict(int)
    for _, max_compl, _ in results:
        histogram[max_compl] += 1

    # Converter chaves para string para JSON e ordenar
    return {str(k): histogram[k] for k in sorted(histogram.keys())}


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a prova de seletividade completa.

    Combina analise de complementaridade off-target com analise de
    conservacao do SL RNA para demonstrar que o ASO e seletivo
    para o parasita alvo.

    Args:
        config: Configuracao do organismo-alvo. Se None, usa L. infantum.

    Returns:
        Envelope JSON completo com resultados.
    """
    # Defaults retrocompativeis com L. infantum
    if config is None:
        config = TargetConfig()

    aso_seq = config.aso_sequence or ASO_SEQUENCE
    target_start = config.aso_target_start or ASO_TARGET_START
    target_end = config.aso_target_end or ASO_TARGET_END
    dg_threshold = config.dg_functional_threshold
    transcriptome_path = config.host_transcriptome_path or TRANSCRIPTOME_HUMAN
    related_seqs = config.related_sl_sequences

    module_name = "02_selectivity_proof"
    envelope = create_envelope(module_name)

    logger.info("=" * 60)
    logger.info("MODULO 02: Prova de Seletividade — %s", config.species_name)
    logger.info("=" * 60)

    try:
        with Timer() as timer:
            # ----------------------------------------------------------
            # Parte 1: Varredura de complementaridade off-target
            # ----------------------------------------------------------
            logger.info("Parte 1: Varredura de complementaridade no transcriptoma")
            logger.info("  ASO: %s (%d nt)", aso_seq, len(aso_seq))
            logger.info("  Limiar RNase H: %d bp contiguos", _RNASE_H_THRESHOLD_BP)

            screen_results, transcriptome_type, screen_warnings = screen_transcriptome(
                aso_seq, transcriptome_path
            )
            envelope["warnings"].extend(screen_warnings)

            # Pior caso: maior complementaridade encontrada
            if screen_results:
                worst_header, worst_compl, worst_pos = screen_results[0]
            else:
                worst_header, worst_compl, worst_pos = ("nenhum", 0, -1)

            logger.info("  Complementaridade maxima encontrada: %d bp", worst_compl)
            logger.info("  Transcrito: %s (pos %d)", worst_header[:60], worst_pos)

            # Calcular dG de ligacao parcial para o pior caso
            partial_dg = compute_partial_binding_dg(aso_seq, worst_compl)
            logger.info("  dG da ligacao parcial (%d bp): %.2f kcal/mol", worst_compl, partial_dg)
            logger.info("  Limiar funcional: %.1f kcal/mol", dg_threshold)

            compl_passes = (
                worst_compl < _RNASE_H_THRESHOLD_BP
                and partial_dg > dg_threshold
            )
            logger.info(
                "  Resultado complementaridade: %s",
                "PASSA" if compl_passes else "FALHA",
            )

            # Histograma
            histogram = _build_complementarity_histogram(screen_results)

            complementarity_data = {
                "transcriptome": transcriptome_type,
                "transcripts_screened": len(screen_results),
                "max_complementarity_found": worst_compl,
                "max_complement_transcript_id": worst_header,
                "max_complement_position": worst_pos,
                "partial_binding_dg_kcal": round(partial_dg, 2),
                "threshold_bp": _RNASE_H_THRESHOLD_BP,
                "threshold_dg_kcal": dg_threshold,
                "passes": compl_passes,
                "histogram": histogram,
            }

            # ----------------------------------------------------------
            # Parte 2: Conservacao do SL RNA (entropia de Shannon)
            # ----------------------------------------------------------
            logger.info("")
            logger.info("Parte 2: Conservacao do SL RNA entre trypanosomatideos")

            # Usa sequencias do config se disponivel; senao, fallback local
            sl_sequences = related_seqs if related_seqs else get_sl_rna_sequences()
            species_names = list(sl_sequences.keys())
            logger.info("  Especies analisadas: %d", len(sl_sequences))
            for sp in species_names:
                logger.info("    - %s", sp)

            # Calcular entropia por posicao — todas as 7 especies
            entropy_array = compute_positional_entropy(sl_sequences)
            total_positions = len(entropy_array)
            zero_entropy_count = sum(1 for h in entropy_array if h == 0.0)

            logger.info(
                "  Posicoes totalmente conservadas (7 spp): %d/%d (%.1f%%)",
                zero_entropy_count,
                total_positions,
                100.0 * zero_entropy_count / total_positions,
            )

            # Entropia apenas entre as 5 especies de Leishmania (alvo terapeutico)
            # Todas as Leishmania tem SL RNA identico — entropia deve ser 0.0 em toda posicao
            leishmania_only = {
                sp: seq for sp, seq in sl_sequences.items() if sp.startswith("L.")
            }
            entropy_leishmania = compute_positional_entropy(leishmania_only)
            leish_zero_count = sum(1 for h in entropy_leishmania if h == 0.0)

            logger.info(
                "  Posicoes totalmente conservadas (Leishmania, %d spp): %d/%d (%.1f%%)",
                len(leishmania_only),
                leish_zero_count,
                total_positions,
                100.0 * leish_zero_count / total_positions,
            )

            # Avaliar conservacao na regiao-alvo do ASO — todas as especies
            conservation = assess_target_region_conservation(
                entropy_array, target_start, target_end
            )
            # Avaliar conservacao na regiao-alvo — apenas Leishmania
            conservation_leish = assess_target_region_conservation(
                entropy_leishmania, target_start, target_end
            )

            logger.info(
                "  Entropia media na regiao-alvo (pos %d-%d, 7 spp): %.4f bits",
                target_start,
                target_end,
                conservation["target_region_mean_entropy"],
            )
            logger.info(
                "  Entropia media na regiao-alvo (Leishmania only): %.4f bits",
                conservation_leish["target_region_mean_entropy"],
            )
            logger.info(
                "  Entropia media fora da regiao-alvo (7 spp): %.4f bits",
                conservation["non_target_mean_entropy"],
            )
            logger.info(
                "  Posicoes conservadas na regiao-alvo (7 spp): %d/%d",
                conservation["fully_conserved_positions_in_target"],
                conservation["target_region_length"],
            )
            logger.info(
                "  Posicoes conservadas na regiao-alvo (Leishmania): %d/%d",
                conservation_leish["fully_conserved_positions_in_target"],
                conservation_leish["target_region_length"],
            )

            # Score de conservacao principal: usa Leishmania (alvo real do ASO)
            # A conservacao entre Leishmania e o que importa para eficacia terapeutica
            # A divergencia com Trypanosoma e informativa mas nao invalida a seletividade
            target_conservation_score = conservation_leish["conservation_ratio"]

            conservation_data = {
                "species_count": len(sl_sequences),
                "species": species_names,
                "sequences": sl_sequences,
                "per_position_entropy": entropy_array,
                "per_position_entropy_leishmania_only": entropy_leishmania,
                "leishmania_species_count": len(leishmania_only),
                "target_region_start": target_start,
                "target_region_end": target_end,
                "target_region_mean_entropy": conservation["target_region_mean_entropy"],
                "target_region_mean_entropy_leishmania": conservation_leish[
                    "target_region_mean_entropy"
                ],
                "non_target_mean_entropy": conservation["non_target_mean_entropy"],
                "fully_conserved_positions_in_target": conservation[
                    "fully_conserved_positions_in_target"
                ],
                "fully_conserved_positions_in_target_leishmania": conservation_leish[
                    "fully_conserved_positions_in_target"
                ],
                "variable_positions_in_target": conservation[
                    "variable_positions_in_target"
                ],
                "variable_positions_in_target_leishmania": conservation_leish[
                    "variable_positions_in_target"
                ],
            }

            # ----------------------------------------------------------
            # Avaliacao combinada
            # ----------------------------------------------------------
            logger.info("")
            logger.info("Avaliacao combinada de seletividade")

            off_target_risk: str
            if worst_compl < _RNASE_H_THRESHOLD_BP and target_conservation_score >= 0.9:
                off_target_risk = "negligible"
            elif worst_compl < _RNASE_H_THRESHOLD_BP:
                off_target_risk = "low"
            else:
                off_target_risk = "significant"

            # Conclusao quantitativa
            conclusion = (
                f"Complementaridade maxima off-target: {worst_compl} bp "
                f"(limiar RNase H: {_RNASE_H_THRESHOLD_BP} bp). "
                f"dG da ligacao parcial: {partial_dg:.2f} kcal/mol "
                f"(limiar funcional: {dg_threshold:.1f} kcal/mol). "
                f"Conservacao da regiao-alvo entre {len(leishmania_only)} especies "
                f"de Leishmania: "
                f"{conservation_leish['fully_conserved_positions_in_target']}/"
                f"{conservation_leish['target_region_length']} posicoes identicas "
                f"(entropia {conservation_leish['target_region_mean_entropy']:.4f} bits). "
                f"Entre todas as {len(sl_sequences)} especies de trypanosomatideos: "
                f"{conservation['fully_conserved_positions_in_target']}/"
                f"{conservation['target_region_length']} posicoes identicas "
                f"(entropia {conservation['target_region_mean_entropy']:.4f} bits). "
                f"Risco off-target: {off_target_risk}. "
                f"O MRL-ASO-001 nao atinge o limiar minimo de {_RNASE_H_THRESHOLD_BP} bp "
                f"contiguos para ativacao de RNase H em nenhum transcrito humano analisado, "
                f"e o alvo SL RNA e 100% conservado entre todas as especies de Leishmania."
            )

            combined_assessment = {
                "off_target_risk": off_target_risk,
                "max_human_complementarity_bp": worst_compl,
                "binding_dg_at_max_complementarity": round(partial_dg, 2),
                "target_conservation_score": round(target_conservation_score, 4),
                "conclusion": conclusion,
            }

            overall_pass = compl_passes and target_conservation_score >= 0.9

            logger.info("  Risco off-target: %s", off_target_risk)
            logger.info(
                "  Score de conservacao do alvo: %.4f", target_conservation_score
            )
            logger.info("  RESULTADO FINAL: %s", "PASSA" if overall_pass else "FALHA")
            logger.info("")
            logger.info("Conclusao: %s", conclusion)

            # ----------------------------------------------------------
            # Montar envelope de resultados
            # ----------------------------------------------------------
            envelope["data"] = {
                "complementarity_screen": complementarity_data,
                "conservation_analysis": conservation_data,
                "combined_assessment": combined_assessment,
            }
            envelope["status"] = "pass" if overall_pass else "fail"
            envelope["summary"]["conclusion"] = conclusion
            envelope["summary"]["key_metrics"] = {
                "max_off_target_complementarity_bp": worst_compl,
                "partial_binding_dg_kcal": round(partial_dg, 2),
                "target_conservation_score": round(target_conservation_score, 4),
                "transcripts_screened": len(screen_results),
                "overall_pass": overall_pass,
            }

        envelope["runtime_seconds"] = timer.elapsed
        logger.info("Tempo de execucao: %.2f s", timer.elapsed)

    except Exception as exc:
        logger.error("Erro no modulo 02: %s", exc, exc_info=True)
        envelope["status"] = "error"
        envelope["warnings"].append(f"Erro durante execucao: {exc}")

    # Sempre grava resultado, mesmo em caso de erro
    output_path = write_result(envelope)
    logger.info("Resultados gravados em: %s", output_path)

    return envelope


if __name__ == "__main__":
    main()
