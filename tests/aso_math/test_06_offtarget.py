"""Testes para o modulo 06 - Off-target screening.

Testa:
  1. Busca de complementaridade encontra matches conhecidos
  2. Scoring classifica risco corretamente
  3. Tratamento de FASTA ausente e de proteinas
  4. Execucao end-to-end do main()
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from aso_math.config import ASO_SEQUENCE, COMPLEMENT

# ---------------------------------------------------------------------------
# Importacoes via importlib (nome do modulo comeca com digito)
# ---------------------------------------------------------------------------

_aligner = importlib.import_module("aso_math.06_offtarget_screen.aligner")
_scorer = importlib.import_module("aso_math.06_offtarget_screen.thermodynamic_scorer")
_run = importlib.import_module("aso_math.06_offtarget_screen.run")

find_max_complementarity = _aligner.find_max_complementarity
find_offtargets = _aligner.find_offtargets
is_nucleotide_fasta = _aligner.is_nucleotide_fasta
parse_fasta = _aligner.parse_fasta
OffTargetHit = _aligner.OffTargetHit

classify_risk = _scorer.classify_risk
score_hits = _scorer.score_hits
ScoredHit = _scorer.ScoredHit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    """Escreve um arquivo FASTA simples."""
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in sequences.items():
            fh.write(">" + header + "\n")
            fh.write(seq + "\n")


def _reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, b) for b in reversed(seq.upper()))


def _complement(seq: str) -> str:
    """Complemento na mesma direcao (nao reverso)."""
    return "".join(COMPLEMENT.get(b, b) for b in seq.upper())


# =========================================================================
# Testes do aligner
# =========================================================================


class TestAligner:
    """Testes para o modulo aligner (busca de complementaridade)."""

    def test_find_max_complementarity_perfect_match(self):
        """Complementaridade perfeita deve retornar comprimento total."""
        aso = "AACCGGTT"
        # Transcript contains the complement (same direction) of the ASO
        transcript = "NNNNN" + _complement(aso) + "NNNNN"
        max_len, pos, region = find_max_complementarity(aso, transcript)
        assert max_len == 8
        assert pos >= 0

    def test_find_max_complementarity_partial_match(self):
        """Match parcial deve retornar comprimento correto."""
        aso = "AACCGGTTAATT"
        partial_compl = _complement("AACCGG")
        transcript = "NNNNNNNN" + partial_compl + "NNNNNNNN"
        max_len, pos, region = find_max_complementarity(aso, transcript)
        assert max_len >= 4

    def test_find_max_complementarity_no_match(self):
        """Sequencias nao-complementares devem retornar 0."""
        # N bases have no complement in the COMPLEMENT dict
        aso = "NNNNNNNNNN"
        transcript = "NNNNNNNNNN"
        max_len, pos, region = find_max_complementarity(aso, transcript)
        assert max_len == 0

    def test_find_offtargets_with_nucleotide_fasta(self, tmp_path):
        """Busca em FASTA de nucleotideos com match conhecido."""
        aso = ASO_SEQUENCE
        rc_18 = _complement(aso[:18])
        fasta_path = tmp_path / "test.fasta"
        _write_fasta(fasta_path, {
            "transcript_with_match": "AAAA" + rc_18 + "AAAA",
            "transcript_no_match": "A" * 100,
        })

        hits = find_offtargets(aso, fasta_path, min_match=15)
        assert len(hits) >= 1
        assert hits[0].match_length >= 15
        assert hits[0].transcript_id == "transcript_with_match"

    def test_find_offtargets_missing_file(self, tmp_path):
        """Arquivo FASTA ausente deve retornar lista vazia."""
        missing = tmp_path / "nonexistent.fasta"
        hits = find_offtargets(ASO_SEQUENCE, missing)
        assert hits == []

    def test_find_offtargets_protein_fasta(self, tmp_path):
        """FASTA de proteinas deve retornar lista vazia."""
        fasta_path = tmp_path / "proteins.fasta"
        _write_fasta(fasta_path, {
            "protein_seq": "MFYWDEKRLHIPQV" * 10,
        })

        hits = find_offtargets(ASO_SEQUENCE, fasta_path)
        assert hits == []

    def test_parse_fasta(self, tmp_path):
        """Parser deve ler headers e sequencias corretamente."""
        fasta_path = tmp_path / "test.fasta"
        _write_fasta(fasta_path, {
            "seq1": "ACGTACGT",
            "seq2": "TTTTGGGG",
        })
        result = parse_fasta(fasta_path)
        assert len(result) == 2
        assert result["seq1"] == "ACGTACGT"
        assert result["seq2"] == "TTTTGGGG"

    def test_is_nucleotide_true(self, tmp_path):
        """FASTA de nucleotideos deve ser detectado."""
        fasta_path = tmp_path / "nuc.fasta"
        _write_fasta(fasta_path, {"s1": "ACGTACGTACGT"})
        assert is_nucleotide_fasta(fasta_path) is True

    def test_is_nucleotide_false(self, tmp_path):
        """FASTA de proteinas deve ser detectado."""
        fasta_path = tmp_path / "prot.fasta"
        _write_fasta(fasta_path, {"s1": "MFYWDEKRLHIPQV"})
        assert is_nucleotide_fasta(fasta_path) is False


# =========================================================================
# Testes do scorer
# =========================================================================


class TestThermodynamicScorer:
    """Testes para o scoring termodinamico e classificacao de risco."""

    def test_classify_risk_safe(self):
        """dG > -15 deve ser classificado como safe."""
        assert classify_risk(-10.0) == "safe"
        assert classify_risk(-14.99) == "safe"
        assert classify_risk(0.0) == "safe"

    def test_classify_risk_monitor(self):
        """dG entre -20 e -15 deve ser classificado como monitor."""
        assert classify_risk(-15.0) == "monitor"
        assert classify_risk(-17.5) == "monitor"
        assert classify_risk(-19.99) == "monitor"

    def test_classify_risk_danger(self):
        """dG <= -20 deve ser classificado como danger."""
        assert classify_risk(-20.0) == "danger"
        assert classify_risk(-25.0) == "danger"
        assert classify_risk(-30.0) == "danger"

    def test_score_hits_ordering(self):
        """Hits devem ser ordenados por dG crescente (mais perigoso primeiro)."""
        hits = [
            OffTargetHit("t1", 10, 0, "A" * 10),
            OffTargetHit("t2", 20, 0, "A" * 20),
        ]
        scored = score_hits(hits, ASO_SEQUENCE)
        assert len(scored) == 2
        assert scored[0].dg_kcal <= scored[1].dg_kcal

    def test_scored_hit_has_risk_level(self):
        """Cada scored hit deve ter um risk_level valido."""
        hits = [OffTargetHit("t1", 12, 0, "A" * 12)]
        scored = score_hits(hits, ASO_SEQUENCE)
        assert scored[0].risk_level in ("safe", "monitor", "danger")

    def test_scored_hit_to_dict(self):
        """to_dict deve retornar dicionario com todas as chaves."""
        hits = [OffTargetHit("t1", 12, 5, "ACGTACGTACGT")]
        scored = score_hits(hits, ASO_SEQUENCE)
        d = scored[0].to_dict()
        expected_keys = {"transcript_id", "match_length", "position",
                         "matched_region", "dg_kcal", "tm_celsius", "risk_level"}
        assert set(d.keys()) == expected_keys


# =========================================================================
# Testes end-to-end
# =========================================================================


class TestModule06EndToEnd:
    """Teste end-to-end do modulo 06."""

    _CACHE: dict = {}

    def _result(self) -> dict:
        if "result" not in self._CACHE:
            self._CACHE["result"] = _run.main()
        return self._CACHE["result"]

    def test_main_runs_without_error(self):
        """main() deve completar sem excecao."""
        result = self._result()
        assert result is not None
        assert isinstance(result, dict)

    def test_envelope_has_required_fields(self):
        """Envelope deve ter todos os campos padrao."""
        result = self._result()
        assert "module" in result
        assert "status" in result
        assert "data" in result
        assert "summary" in result
        assert result["module"] == "06_offtarget_screen"

    def test_status_not_error(self):
        """Status nao deve ser error."""
        result = self._result()
        assert result["status"] != "error"

    def test_data_has_risk_summary(self):
        """Dados devem conter resumo de risco."""
        result = self._result()
        data = result["data"]
        assert "risk_summary" in data
        assert "safe" in data["risk_summary"]
        assert "monitor" in data["risk_summary"]
        assert "danger" in data["risk_summary"]

    def test_data_has_transcripts_screened(self):
        """Dados devem reportar quantos transcritos foram varridos."""
        result = self._result()
        assert result["data"]["transcripts_screened"] > 0

    def test_null_model_fallback(self):
        """Sem transcriptoma canino, deve usar modelo nulo."""
        result = self._result()
        assert result["data"]["transcriptome_source"] == "null_model"

    def test_summary_has_key_metrics(self):
        """Summary deve conter key_metrics com campos esperados."""
        result = self._result()
        km = result["summary"]["key_metrics"]
        assert "transcripts_screened" in km
        assert "total_hits" in km
        assert "danger_hits" in km
        assert "overall_safe" in km
