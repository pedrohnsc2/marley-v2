"""Testes unitarios para aso_math/organisms.py.

Valida o registro de organismos e a construcao de TargetConfig:
- list_organisms() retorna todos os organismos registrados
- get_target_config() constroi configs corretas por slug
- Erro para organismo desconhecido
- Override de parametros
"""

import pytest

from aso_math.organisms import get_target_config, list_organisms


# ---------------------------------------------------------------------------
# 1. Listagem de organismos
# ---------------------------------------------------------------------------


def test_list_organisms_not_empty():
    """Registro deve conter pelo menos 15 organismos."""
    orgs = list_organisms()
    assert len(orgs) >= 15


# ---------------------------------------------------------------------------
# 2. Leishmania infantum (organismo principal)
# ---------------------------------------------------------------------------


def test_get_leishmania_infantum():
    """Config de L. infantum deve ter valores conhecidos do pipeline."""
    config = get_target_config("leishmania_infantum")
    assert config.species_name == "Leishmania infantum"
    assert config.sl_length == 39
    assert config.known_tm == 68.48
    assert config.known_dg == -27.97


# ---------------------------------------------------------------------------
# 3. Trypanosoma cruzi
# ---------------------------------------------------------------------------


def test_get_trypanosoma_cruzi():
    """T. cruzi deve ter SL de 39 nt e >= 10 sequencias relacionadas."""
    config = get_target_config("trypanosoma_cruzi")
    assert config.species_name == "Trypanosoma cruzi"
    assert config.sl_length == 39
    assert len(config.related_sl_sequences) >= 10


# ---------------------------------------------------------------------------
# 4. Nematodeo (Brugia malayi)
# ---------------------------------------------------------------------------


def test_get_nematode():
    """B. malayi deve ter SL1 de 22 nt e grupo taxonomico nematoda."""
    config = get_target_config("brugia_malayi")
    assert config.sl_length == 22
    assert config.taxonomic_group == "nematoda"


# ---------------------------------------------------------------------------
# 5. Organismo desconhecido
# ---------------------------------------------------------------------------


def test_unknown_organism_raises():
    """Slug invalido deve levantar ValueError com mensagem em portugues."""
    with pytest.raises(ValueError, match="Organismo desconhecido"):
        get_target_config("fake_species")


# ---------------------------------------------------------------------------
# 6. Override de parametros
# ---------------------------------------------------------------------------


def test_override_parameter():
    """Override via kwargs deve alterar apenas o campo especificado."""
    config = get_target_config("leishmania_infantum", lna_5prime=2)
    assert config.lna_5prime == 2
    # Demais campos devem permanecer inalterados
    assert config.species_name == "Leishmania infantum"
