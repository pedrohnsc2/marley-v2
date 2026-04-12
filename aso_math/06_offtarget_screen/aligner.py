"""Alinhamento de complementaridade off-target por janela deslizante.

Busca regioes de complementaridade contigua entre o ASO e transcritos
de um arquivo FASTA, reutilizando a logica do modulo 02 (find_max_complementarity).

Sem dependencias externas — parser FASTA e busca implementados em Python puro.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from aso_math.config import COMPLEMENT


# ---------------------------------------------------------------------------
# Tipos de retorno
# ---------------------------------------------------------------------------


@dataclass
class OffTargetHit:
    """Um hit de complementaridade off-target."""

    transcript_id: str
    match_length: int
    position: int
    matched_region: str

    def to_dict(self) -> dict[str, Any]:
        """Serializa para dicionario JSON-compativel."""
        return {
            "transcript_id": self.transcript_id,
            "match_length": self.match_length,
            "position": self.position,
            "matched_region": self.matched_region,
        }


# ---------------------------------------------------------------------------
# Deteccao de tipo FASTA (nucleotideo vs proteina)
# ---------------------------------------------------------------------------


def is_nucleotide_fasta(fasta_path: Path, sample_lines: int = 20) -> bool:
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
                if protein_only_chars & set(stripped):
                    return False
                lines_checked += 1
                if lines_checked >= sample_lines:
                    break
    except OSError:
        return False

    return lines_checked > 0


# ---------------------------------------------------------------------------
# Parser FASTA
# ---------------------------------------------------------------------------


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
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq).upper()
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)

    if current_header is not None:
        sequences[current_header] = "".join(current_seq).upper()

    return sequences


# ---------------------------------------------------------------------------
# Busca de complementaridade (janela deslizante)
# ---------------------------------------------------------------------------


def _reverse_complement(seq: str) -> str:
    """Complemento reverso de uma sequencia DNA/RNA."""
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq.upper()))


def find_max_complementarity(
    aso_seq: str, transcript_seq: str
) -> tuple[int, int, str]:
    """Encontra o comprimento maximo de complementaridade contigua.

    Desliza o ASO (e seu complemento reverso) ao longo de um transcrito.
    Em cada posicao, conta o maior run de bases complementares contiguas.

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        transcript_seq: Sequencia do transcrito.

    Returns:
        Tupla (max_length, best_position, matched_region) onde position e
        a posicao no transcrito onde o melhor alinhamento comeca.
    """
    aso_upper = aso_seq.upper()
    transcript_upper = transcript_seq.upper()
    aso_rc = _reverse_complement(aso_upper)

    global_max_length = 0
    global_best_pos = -1

    for probe in (aso_upper, aso_rc):
        probe_len = len(probe)
        transcript_len = len(transcript_upper)

        for offset in range(-(probe_len - 1), transcript_len):
            t_start = max(0, offset)
            p_start = max(0, -offset)
            overlap = min(probe_len - p_start, transcript_len - t_start)

            if overlap < 1:
                continue

            current_run = 0
            best_run = 0
            best_run_pos = t_start

            for i in range(overlap):
                t_base = transcript_upper[t_start + i]
                p_base = probe[p_start + i]

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

    matched_region = ""
    if global_max_length > 0 and global_best_pos >= 0:
        matched_region = transcript_upper[
            global_best_pos : global_best_pos + global_max_length
        ]

    return global_max_length, global_best_pos, matched_region


# ---------------------------------------------------------------------------
# API publica
# ---------------------------------------------------------------------------


def find_offtargets(
    aso_seq: str,
    fasta_path: Path,
    min_match: int = 15,
) -> list[OffTargetHit]:
    """Busca hits off-target em um transcriptoma FASTA.

    Para cada transcrito no FASTA, computa a complementaridade maxima
    contigua com o ASO. Retorna apenas hits com complementaridade
    >= min_match bases.

    Detecta automaticamente se o FASTA contem proteinas (e pula nesse caso).
    Retorna lista vazia se o arquivo nao existir ou for proteina.

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        fasta_path: Caminho para arquivo FASTA do transcriptoma.
        min_match: Comprimento minimo de complementaridade para reportar
                   como hit (default: 15 bp).

    Returns:
        Lista de OffTargetHit ordenada por match_length decrescente.
    """
    if not fasta_path.exists():
        return []

    if not is_nucleotide_fasta(fasta_path):
        return []

    sequences = parse_fasta(fasta_path)
    hits: list[OffTargetHit] = []

    for header, seq in sequences.items():
        max_compl, position, matched_region = find_max_complementarity(
            aso_seq, seq
        )
        if max_compl >= min_match:
            hits.append(
                OffTargetHit(
                    transcript_id=header,
                    match_length=max_compl,
                    position=position,
                    matched_region=matched_region,
                )
            )

    hits.sort(key=lambda h: h.match_length, reverse=True)
    return hits
