"""Registro de organismos com SL RNA para o pipeline ASO.

Cada entrada contem a sequencia do SL RNA, parametros biologicos
e sequencias de especies relacionadas para analise de conservacao.

Uso:
    from aso_math.organisms import get_target_config, list_organisms
    config = get_target_config("trypanosoma_cruzi")
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from aso_math.target_config import TargetConfig


# ---------------------------------------------------------------------------
# Sequencias SL RNA compartilhadas entre kinetoplastideos
# Ref: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
# Ref: Milhausen M et al. (1984) Cell 38(3):721-729
# Ref: Sturm NR et al. (1999) Mol Biochem Parasitol 104(1):69-80
# ---------------------------------------------------------------------------

_KINETOPLASTID_SL_SEQUENCES: dict[str, str] = {
    "Leishmania infantum": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Leishmania donovani": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Leishmania major": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Leishmania braziliensis": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Leishmania mexicana": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Leishmania amazonensis": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
    "Trypanosoma brucei": "AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG",
    "Trypanosoma cruzi": "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG",
    "Trypanosoma vivax": "AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG",
    "Crithidia fasciculata": "AACTAACGCTATATAAGTATCAGTTTCTGTACTATATCG",
    "Leptomonas seymouri": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
}

# ---------------------------------------------------------------------------
# Sequencias SL1 de nematodeos
# Ref: Bektesh S et al. (1988) Genes Dev 2:1277-1283
# O SL1 de 22 nt e conservado em todo o filo Nematoda
# ---------------------------------------------------------------------------

_NEMATODE_SL1_SEQUENCES: dict[str, str] = {
    "Caenorhabditis elegans": "GGTTTAATTACCCAAGTTTGAG",
    "Caenorhabditis briggsae": "GGTTTAATTACCCAAGTTTGAG",
    "Ascaris lumbricoides": "GGTTTAATTACCCAAGTTTGAG",
    "Brugia malayi": "GGTTTAATTACCCAAGTTTGAG",
    "Strongyloides stercoralis": "GGTTTAATTACCCAAGTTTGAG",
    "Necator americanus": "GGTTTAATTACCCAAGTTTGAG",
    "Onchocerca volvulus": "GGTTTAATTACCCAAGTTTGAG",
    "Wuchereria bancrofti": "GGTTTAATTACCCAAGTTTGAG",
    "Dirofilaria immitis": "GGTTTAATTACCCAAGTTTGAG",
}

# ---------------------------------------------------------------------------
# Registro de organismos
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, dict] = {
    # =====================================================================
    # KINETOPLASTIDA
    # =====================================================================
    "leishmania_infantum": {
        "species_name": "Leishmania infantum",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "leishmaniose visceral",
        "sl_sequence": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "aso_sequence": "ACAGAAACTGATACTTATATAGCGT",
        "aso_name": "MRL-ASO-001",
        "aso_target_start": 5,
        "aso_target_end": 30,
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 12.0,
        "sl_copy_number": 150,
        "divergence_time_mya": 350.0,
        "known_tm": 68.48,
        "known_dg": -27.97,
        "known_gc": 0.32,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    # Ref: Fernandez-Moya SM et al. (2009) Parasitology 136(8):851-861
    "leishmania_braziliensis": {
        "species_name": "Leishmania braziliensis",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "leishmaniose mucocutanea",
        "sl_sequence": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 12.0,
        "sl_copy_number": 150,
        "divergence_time_mya": 350.0,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    # Ref: Zilka A et al. (2001) J Biol Chem 276(50):47261-47267
    "leishmania_donovani": {
        "species_name": "Leishmania donovani",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "leishmaniose visceral (calazar)",
        "sl_sequence": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 12.0,
        "sl_copy_number": 150,
        "divergence_time_mya": 350.0,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    # Ref: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
    "leishmania_major": {
        "species_name": "Leishmania major",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "leishmaniose cutanea",
        "sl_sequence": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 12.0,
        "sl_copy_number": 150,
        "divergence_time_mya": 350.0,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    "trypanosoma_brucei": {
        "species_name": "Trypanosoma brucei",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "doenca do sono",
        "sl_sequence": "AACTAACGCTATTATTAGAACAGTTTCTGTACTATATTG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 6.0,
        "sl_copy_number": 200,
        "divergence_time_mya": 350.0,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    "trypanosoma_cruzi": {
        "species_name": "Trypanosoma cruzi",
        "taxonomic_group": "kinetoplastida",
        "disease_name": "doenca de Chagas",
        "sl_sequence": "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 24.0,
        "sl_copy_number": 200,
        "divergence_time_mya": 350.0,
        "related_sl_sequences": _KINETOPLASTID_SL_SEQUENCES,
    },
    # =====================================================================
    # NEMATODA
    # Ref SL1: Bektesh S et al. (1988) Genes Dev 2:1277-1283
    # =====================================================================
    "ascaris_lumbricoides": {
        "species_name": "Ascaris lumbricoides",
        "taxonomic_group": "nematoda",
        "disease_name": "ascaridiase",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 1440.0,  # ~60 dias
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    "brugia_malayi": {
        "species_name": "Brugia malayi",
        "taxonomic_group": "nematoda",
        "disease_name": "filariose linfatica",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 720.0,  # ~30 dias (ciclo do verme adulto)
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    "dirofilaria_immitis": {
        "species_name": "Dirofilaria immitis",
        "taxonomic_group": "nematoda",
        "disease_name": "dirofilariose (verme do coracao canino)",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 4380.0,  # ~6 meses
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    # Ref: Blaxter ML et al. (1998) Nature 392:71-75
    "necator_americanus": {
        "species_name": "Necator americanus",
        "taxonomic_group": "nematoda",
        "disease_name": "ancilostomiase",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 1440.0,  # ~60 dias
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    "onchocerca_volvulus": {
        "species_name": "Onchocerca volvulus",
        "taxonomic_group": "nematoda",
        "disease_name": "oncocercose (cegueira dos rios)",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 8760.0,  # ~1 ano (verme adulto vive 10-14 anos)
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    # Ref: Guiliano DB & Bhopale MK (2003) Mol Biochem Parasitol 130:93-103
    "strongyloides_stercoralis": {
        "species_name": "Strongyloides stercoralis",
        "taxonomic_group": "nematoda",
        "disease_name": "estrongiloidiase",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 504.0,  # ~21 dias
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    # Ref: Blaxter ML et al. (1998) Nature 392:71-75
    "wuchereria_bancrofti": {
        "species_name": "Wuchereria bancrofti",
        "taxonomic_group": "nematoda",
        "disease_name": "filariose linfatica (elefantiase)",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 8760.0,  # ~1 ano
        "sl_copy_number": 110,
        "length_scan_min": 14,
        "length_scan_max": 22,
        "divergence_time_mya": 500.0,
        "neutral_rate_per_site_per_year": 2.0e-9,
        "related_sl_sequences": _NEMATODE_SL1_SEQUENCES,
    },
    # =====================================================================
    # PLATYHELMINTHES
    # Ref SL: Davis RE et al. (1995) J Biol Chem 270(37):21813-21819
    # =====================================================================
    "schistosoma_haematobium": {
        "species_name": "Schistosoma haematobium",
        "taxonomic_group": "platyhelminthes",
        "disease_name": "esquistossomose urogenital",
        "sl_sequence": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",  # 36 nt
        "mutation_rate": 3.0e-9,
        "generation_time_hours": 720.0,  # ~30 dias
        "sl_copy_number": 100,
        "length_scan_min": 16,
        "length_scan_max": 27,
        "divergence_time_mya": 300.0,
        "neutral_rate_per_site_per_year": 3.0e-9,
        "related_sl_sequences": {
            "Schistosoma mansoni": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",
            "Schistosoma haematobium": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",
        },
    },
    # Ref: Davis RE et al. (1995) J Biol Chem 270(37):21813-21819
    "schistosoma_mansoni": {
        "species_name": "Schistosoma mansoni",
        "taxonomic_group": "platyhelminthes",
        "disease_name": "esquistossomose intestinal",
        "sl_sequence": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",  # 36 nt
        "mutation_rate": 3.0e-9,
        "generation_time_hours": 720.0,  # ~30 dias
        "sl_copy_number": 100,
        "length_scan_min": 16,
        "length_scan_max": 27,
        "divergence_time_mya": 300.0,
        "neutral_rate_per_site_per_year": 3.0e-9,
        "related_sl_sequences": {
            "Schistosoma mansoni": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",
            "Schistosoma japonicum": "AACCGTCACGGTTTTACTCTTGTGATTTGTTGCATG",
        },
    },
}


def get_target_config(slug: str, **overrides: Any) -> TargetConfig:
    """Constroi um TargetConfig a partir do registro.

    Args:
        slug: Identificador do organismo (ex: "trypanosoma_cruzi").
        **overrides: Campos a sobrescrever (ex: aso_sequence="...", lna_5prime=3).

    Returns:
        TargetConfig configurado para o organismo.

    Raises:
        ValueError: Se o organismo nao estiver no registro.
    """
    if slug not in _REGISTRY:
        available = ", ".join(sorted(_REGISTRY.keys()))
        raise ValueError(
            f"Organismo desconhecido: '{slug}'. Disponiveis: {available}"
        )
    params: dict[str, Any] = {"organism_slug": slug, **_REGISTRY[slug]}
    params.update(overrides)
    return TargetConfig(**params)


def list_organisms() -> list[dict[str, str]]:
    """Retorna lista de organismos disponiveis no registro."""
    result = []
    for slug, data in _REGISTRY.items():
        result.append({
            "slug": slug,
            "species": data["species_name"],
            "disease": data.get("disease_name", ""),
            "group": data.get("taxonomic_group", ""),
            "sl_length": str(len(data["sl_sequence"])),
        })
    return result
