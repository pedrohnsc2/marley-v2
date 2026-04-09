"""Testes unitarios para aso_math/thermo.py.

Valida funcoes termodinamicas do modelo nearest-neighbor (SantaLucia 1998):
- reverse_complement, gc_content
- compute_tm, compute_dg
- compute_hairpin_dg, compute_self_dimer_dg

Valores de referencia: MRL-ASO-001 (Tm=68.48, dG=-27.97).
"""

import pytest

from aso_math.thermo import (
    compute_dg,
    compute_hairpin_dg,
    compute_self_dimer_dg,
    compute_tm,
    gc_content,
    reverse_complement,
)


# ---------------------------------------------------------------------------
# 1. Valores conhecidos — MRL-ASO-001
# Referencia: pipeline original (aso_candidates.json, config.py)
# ---------------------------------------------------------------------------


def test_compute_tm_mrl_aso_001():
    """Tm do MRL-ASO-001 deve ser 68.48 C (valor validado no pipeline)."""
    assert compute_tm("ACAGAAACTGATACTTATATAGCGT") == 68.48


def test_compute_dg_mrl_aso_001():
    """dG do MRL-ASO-001 deve ser -27.97 kcal/mol (valor validado no pipeline)."""
    assert compute_dg("ACAGAAACTGATACTTATATAGCGT") == -27.97


# ---------------------------------------------------------------------------
# 2. Complemento reverso
# ---------------------------------------------------------------------------


def test_reverse_complement_basic():
    """Complemento reverso basico: ATGC -> GCAT."""
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_with_u():
    """U (RNA) deve ser tratado como T para complemento: AUGC -> GCAT."""
    assert reverse_complement("AUGC") == "GCAT"


def test_reverse_complement_palindrome():
    """Sequencia palindromica (sitio EcoRI) deve retornar ela mesma."""
    seq = "GAATTC"
    assert reverse_complement(seq) == "GAATTC"


# ---------------------------------------------------------------------------
# 3. Conteudo GC
# ---------------------------------------------------------------------------


def test_gc_content_all_gc():
    """Sequencia 100% GC deve retornar 1.0."""
    assert gc_content("GGCC") == 1.0


def test_gc_content_all_at():
    """Sequencia 0% GC (so AT) deve retornar 0.0."""
    assert gc_content("AATT") == 0.0


def test_gc_content_mrl_aso_001():
    """GC do MRL-ASO-001 deve ser ~0.32 (8 GC de 25 nt)."""
    assert abs(gc_content("ACAGAAACTGATACTTATATAGCGT") - 0.32) < 0.01


def test_gc_content_empty():
    """Sequencia vazia deve retornar 0.0 (sem divisao por zero)."""
    assert gc_content("") == 0.0


# ---------------------------------------------------------------------------
# 4. Tm — casos limite
# ---------------------------------------------------------------------------


def test_compute_tm_short_sequence_wallace_rule():
    """Para seq < 14 nt, usa regra de Wallace: Tm = 2*(A+T) + 4*(G+C)."""
    seq = "ATGCATGCATGC"  # 12 nt, 6 AT + 6 GC
    expected = 2 * 6 + 4 * 6  # = 36
    assert compute_tm(seq) == float(expected)


def test_compute_tm_long_gc_rich():
    """Sequencia rica em GC deve ter Tm maior que sequencia rica em AT."""
    gc_rich = "GCGCGCGCGCGCGCGCGC"  # 18 nt, 100% GC
    at_rich = "ATATATATATATATATAT"  # 18 nt, 0% GC
    assert compute_tm(gc_rich) > compute_tm(at_rich)


# ---------------------------------------------------------------------------
# 5. dG — propriedades
# ---------------------------------------------------------------------------


def test_compute_dg_negative():
    """dG de hibridizacao deve ser negativo (ligacao favoravel)."""
    assert compute_dg("ACAGAAACTGATACTTATATAGCGT") < 0


def test_compute_dg_longer_more_negative():
    """Sequencia mais longa deve ter dG mais negativo (ligacao mais forte)."""
    short = "ACAGAAACTGATACT"  # 15 nt
    long = "ACAGAAACTGATACTTATATAGCGT"  # 25 nt
    assert compute_dg(long) < compute_dg(short)


# ---------------------------------------------------------------------------
# 6. Hairpin dG
# ---------------------------------------------------------------------------


def test_hairpin_dg_no_hairpin():
    """MRL-ASO-001 nao deve formar hairpin forte (dG >= -5.0)."""
    assert compute_hairpin_dg("ACAGAAACTGATACTTATATAGCGT") >= -5.0


def test_hairpin_dg_palindrome():
    """Sequencia palindromica deve retornar float valido."""
    seq = "AATTAATTAATTAATT"
    dg = compute_hairpin_dg(seq)
    assert isinstance(dg, float)


# ---------------------------------------------------------------------------
# 7. Self-dimer dG
# ---------------------------------------------------------------------------


def test_self_dimer_dg_returns_float():
    """Self-dimer dG deve ser float <= 0.0."""
    dg = compute_self_dimer_dg("ACAGAAACTGATACTTATATAGCGT")
    assert isinstance(dg, float)
    assert dg <= 0.0  # sempre <= 0 (0.0 = sem dimerizacao)
