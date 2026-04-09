"""Testes unitarios para aso_math/target_config.py.

Valida o dataclass TargetConfig:
- Valores default (L. infantum)
- Auto-calculo de sl_length via __post_init__
- Imutabilidade (frozen=True)
- Configuracao customizada
"""

import pytest

from aso_math.target_config import TargetConfig


# ---------------------------------------------------------------------------
# 1. Configuracao default (L. infantum)
# ---------------------------------------------------------------------------


def test_default_config():
    """Defaults devem corresponder a L. infantum."""
    config = TargetConfig()
    assert config.organism_slug == "leishmania_infantum"
    assert config.sl_length == 39
    assert len(config.sl_sequence) == 39


# ---------------------------------------------------------------------------
# 2. Auto-calculo de sl_length
# ---------------------------------------------------------------------------


def test_auto_sl_length():
    """sl_length deve ser auto-calculado a partir de sl_sequence."""
    config = TargetConfig(sl_sequence="ATGCATGC")
    assert config.sl_length == 8


# ---------------------------------------------------------------------------
# 3. Imutabilidade (frozen dataclass)
# ---------------------------------------------------------------------------


def test_frozen():
    """TargetConfig e frozen — atribuicao deve levantar excecao."""
    config = TargetConfig()
    with pytest.raises((AttributeError,)):
        config.organism_slug = "other"  # type: ignore[misc]


# ---------------------------------------------------------------------------
# 4. Configuracao customizada
# ---------------------------------------------------------------------------


def test_custom_organism():
    """Construcao com parametros customizados deve funcionar."""
    config = TargetConfig(
        organism_slug="test",
        species_name="Test species",
        sl_sequence="GGTTTAATTACCCAAGTTTGAG",
    )
    assert config.sl_length == 22
    assert config.organism_slug == "test"
    assert config.species_name == "Test species"
