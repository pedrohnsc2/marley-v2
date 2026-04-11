"""Execucao do modulo 05_rosettafold — Analise estrutural do complexo ASO:SL RNA.

Pipeline computacional que substitui RoseTTAFold2NA por analise geometrica
direta. Constroi modelos 3D, analisa geometria e decompoe energia de ligacao.

Uso:
    python -m marley_ai.05_rosettafold.run
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import importlib as _il

from marley_ai.config import AIModuleConfig, detect_device
from marley_ai.envelope import Timer, create_envelope, write_result
from marley_ai.registry import register

# Importacao via importlib — "05" nao e valido como identificador Python
_analyzer = _il.import_module("marley_ai.05_rosettafold.analyzer")
_builder = _il.import_module("marley_ai.05_rosettafold.builder")
_cfg = _il.import_module("marley_ai.05_rosettafold.config")
_scoring = _il.import_module("marley_ai.05_rosettafold.scoring")

StructureAnalysisResult = _analyzer.StructureAnalysisResult
analyze_structure = _analyzer.analyze_structure
parse_pdb = _analyzer.parse_pdb

atoms_to_pdb = _builder.atoms_to_pdb
build_aso_sl_complex = _builder.build_aso_sl_complex
save_pdb = _builder.save_pdb

ASO_PDB = _cfg.ASO_PDB
ASO_SEQUENCE = _cfg.ASO_SEQUENCE
ASO_TARGET_END = _cfg.ASO_TARGET_END
ASO_TARGET_START = _cfg.ASO_TARGET_START
LNA_POSITIONS = _cfg.LNA_POSITIONS
OUTPUT_RESULTS_DIR = _cfg.OUTPUT_RESULTS_DIR
OUTPUT_STRUCTURES_DIR = _cfg.OUTPUT_STRUCTURES_DIR
SL_RNA_PDB = _cfg.SL_RNA_PDB
SL_RNA_SEQUENCE = _cfg.SL_RNA_SEQUENCE
SL_SEQUENCE = _cfg.SL_SEQUENCE
StructuralAnalysisConfig = _cfg.StructuralAnalysisConfig

BindingEnergyResult = _scoring.BindingEnergyResult
decompose_binding_energy = _scoring.decompose_binding_energy


# ---------------------------------------------------------------------------
# Registro do modulo no sistema de registro global
# ---------------------------------------------------------------------------

@register("05_rosettafold")
class StructuralAnalysisModule:
    """Modulo de analise estrutural do complexo ASO:SL RNA.

    Substitui a dependencia de RoseTTAFold2NA por um pipeline computacional
    que constroi modelos geometricos, analisa estruturas existentes, e
    decompoe a energia de ligacao em componentes fisicos.
    """

    def __init__(self) -> None:
        self._config: StructuralAnalysisConfig | None = None
        self._result: dict[str, Any] | None = None

    def configure(self, config: Any) -> None:
        """Recebe configuracao e prepara o modulo."""
        if isinstance(config, StructuralAnalysisConfig):
            self._config = config
        else:
            self._config = StructuralAnalysisConfig()

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica dependencias: sequencias e (opcionalmente) PDBs existentes."""
        missing: list[str] = []

        # Sequencias sao hardcoded — sempre disponiveis
        if not ASO_SEQUENCE:
            missing.append("ASO_SEQUENCE")
        if not SL_SEQUENCE:
            missing.append("SL_SEQUENCE")

        # PDBs existentes sao opcionais (analise comparativa)
        existing_pdbs: list[str] = []
        if ASO_PDB.exists():
            existing_pdbs.append(str(ASO_PDB))
        if SL_RNA_PDB.exists():
            existing_pdbs.append(str(SL_RNA_PDB))

        return {
            "valid": len(missing) == 0,
            "missing": missing,
            "existing_pdbs": existing_pdbs,
        }

    def run(self) -> dict[str, Any]:
        """Executa o pipeline completo de analise estrutural."""
        return main(self._config)

    def get_dependencies(self) -> list[str]:
        """Dependencias: nenhuma obrigatoria. aso_math e opcional."""
        return []


# ---------------------------------------------------------------------------
# Funcoes auxiliares para formatacao de resultados
# ---------------------------------------------------------------------------

def _analysis_to_dict(analysis: StructureAnalysisResult) -> dict[str, Any]:
    """Converte StructureAnalysisResult em dicionario serializavel."""
    return {
        "helical_parameters": {
            "mean_rise_angstrom": round(analysis.helical_params.mean_rise, 2),
            "std_rise": round(analysis.helical_params.std_rise, 2),
            "mean_twist_degrees": round(analysis.helical_params.mean_twist, 1),
            "std_twist": round(analysis.helical_params.std_twist, 1),
            "bp_per_turn": round(analysis.helical_params.bp_per_turn, 1),
            "n_base_pairs": analysis.helical_params.n_base_pairs,
            "helix_length_angstrom": round(analysis.helical_params.helix_length, 1),
            "helix_form": analysis.helical_params.helix_form,
        },
        "groove_analysis": {
            "major_groove_width_angstrom": round(
                analysis.groove_analysis.major_groove_width, 1
            ),
            "minor_groove_width_angstrom": round(
                analysis.groove_analysis.minor_groove_width, 1
            ),
            "minor_to_major_ratio": round(
                analysis.groove_analysis.minor_to_major_ratio, 2
            ),
        },
        "hydrogen_bonds": {
            "total": analysis.n_total_hbonds,
            "watson_crick": analysis.n_watson_crick,
            "wobble": analysis.n_wobble,
            "non_canonical": (
                analysis.n_total_hbonds - analysis.n_watson_crick - analysis.n_wobble
            ),
        },
        "rmsd": {
            "vs_existing_pdb_angstrom": round(analysis.rmsd_vs_existing, 2),
        },
        "sasa": {
            "free_rna_angstrom2": round(analysis.sasa_free, 1),
            "bound_complex_angstrom2": round(analysis.sasa_bound, 1),
            "delta_sasa_angstrom2": round(analysis.delta_sasa, 1),
        },
        "rnase_h_accessibility": {
            "accessible": analysis.rnase_h_accessible,
            "score": round(analysis.rnase_h_score, 2),
        },
        "structural_quality": analysis.structural_quality,
    }


def _energy_to_dict(energy: BindingEnergyResult) -> dict[str, Any]:
    """Converte BindingEnergyResult em dicionario serializavel."""
    components = {
        "nearest_neighbor": {
            "dg_kcal_mol": round(energy.nn_thermodynamic.value_kcal, 2),
            "description": energy.nn_thermodynamic.description,
        },
        "base_stacking": {
            "dg_kcal_mol": round(energy.stacking_energy.value_kcal, 2),
            "n_contributions": energy.stacking_energy.n_contributions,
        },
        "hydrogen_bonding": {
            "dg_kcal_mol": round(energy.hbond_energy.value_kcal, 2),
            "n_contributions": energy.hbond_energy.n_contributions,
        },
        "electrostatic": {
            "dg_kcal_mol": round(energy.electrostatic_energy.value_kcal, 2),
            "description": energy.electrostatic_energy.description,
        },
        "desolvation_penalty": {
            "dg_kcal_mol": round(energy.solvation_penalty.value_kcal, 2),
        },
        "lna_stabilization": {
            "dg_kcal_mol": round(energy.lna_contribution.value_kcal, 2),
            "n_positions": energy.lna_contribution.n_contributions,
            "description": energy.lna_contribution.description,
        },
        "ps_backbone": {
            "dg_kcal_mol": round(energy.ps_contribution.value_kcal, 2),
            "n_positions": energy.ps_contribution.n_contributions,
        },
    }

    return {
        "components": components,
        "totals": {
            "dg_total_kcal_mol": round(energy.total_dg, 2),
            "dh_total_kcal_mol": round(energy.total_dh, 2),
            "ds_total_cal_mol_k": round(energy.total_ds, 2),
            "predicted_tm_celsius": round(energy.predicted_tm, 1),
        },
        "experimental_comparison": {
            "experimental_dg_kcal_mol": energy.experimental_dg,
            "experimental_tm_celsius": energy.experimental_tm,
            "dg_deviation_percent": round(energy.dg_deviation, 1),
            "tm_deviation_celsius": round(energy.tm_deviation, 1),
        },
        "classification": {
            "is_functional": energy.is_functional,
            "binding_class": energy.binding_class,
            "threshold_kcal_mol": -15.0,
        },
        "lna_vs_dna": {
            k: round(v, 3) if isinstance(v, float) else v
            for k, v in energy.lna_vs_dna_comparison.items()
        },
    }


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------

def main(config: StructuralAnalysisConfig | None = None) -> dict[str, Any]:
    """Executa pipeline completo de analise estrutural do complexo ASO:SL RNA.

    Etapas:
    1. Construir modelos 3D (SL RNA isolado + duplex + complexo)
    2. Analisar geometria (helice, sulcos, H-bonds, SASA)
    3. Decompor energia de ligacao (NN, stacking, H-bond, eletrostatica)
    4. Comparar com PDBs existentes (se disponiveis)
    5. Salvar PDB + resultados JSON

    Args:
        config: Configuracao do modulo. Se None, usa defaults.

    Returns:
        Envelope com resultados completos da analise.
    """
    if config is None:
        config = StructuralAnalysisConfig()

    envelope = create_envelope("05_rosettafold", version="1.0.0")
    envelope["device"] = detect_device()
    envelope["dependencies"] = []

    warnings: list[str] = []

    with Timer() as timer:
        # ==================================================================
        # ETAPA 1: Construcao de modelos 3D
        # ==================================================================
        print("[05_rosettafold] Etapa 1: Construindo modelos 3D...")

        models = build_aso_sl_complex(config)
        sl_rna_free = models["sl_rna_free"]
        duplex = models["duplex"]
        complex_atoms = models["complex"]

        print(f"  SL RNA isolado: {len(sl_rna_free)} atomos")
        print(f"  Duplex ASO:RNA: {len(duplex)} atomos")
        print(f"  Complexo completo: {len(complex_atoms)} atomos")

        # Salvar PDBs
        OUTPUT_STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)

        pdb_sl_free = save_pdb(
            sl_rna_free,
            OUTPUT_STRUCTURES_DIR / "sl_rna_free_model.pdb",
            title="SL RNA L. infantum - Free Form (39 nt)",
            remarks=[
                f"Sequence: {SL_RNA_SEQUENCE}",
                f"Length: {len(SL_RNA_SEQUENCE)} nt",
                "Idealized A-form single strand helix",
            ],
        )

        pdb_duplex = save_pdb(
            duplex,
            OUTPUT_STRUCTURES_DIR / "aso_sl_duplex_model.pdb",
            title="MRL-ASO-001 : SL RNA Duplex (24 bp)",
            remarks=[
                f"RNA (chain A): {SL_RNA_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END]}",
                f"ASO (chain B): {ASO_SEQUENCE}",
                f"LNA positions: {LNA_POSITIONS}",
                "A-form RNA:DNA hybrid with PS backbone",
                "B-factor coding: 20=LNA, 10=PS, 0=unmodified",
            ],
        )

        pdb_complex = save_pdb(
            complex_atoms,
            OUTPUT_STRUCTURES_DIR / "aso_sl_complex_model.pdb",
            title="MRL-ASO-001 : SL RNA Full Complex",
            remarks=[
                f"SL RNA (chain A): {SL_RNA_SEQUENCE} (39 nt)",
                f"ASO (chain B): {ASO_SEQUENCE} (24 nt)",
                f"Hybridization region: positions {ASO_TARGET_START}-{ASO_TARGET_END}",
                "5' flank + duplex + 3' flank",
            ],
        )

        artifacts = [
            str(pdb_sl_free),
            str(pdb_duplex),
            str(pdb_complex),
        ]

        print(f"  PDBs salvos em {OUTPUT_STRUCTURES_DIR}/")

        # ==================================================================
        # ETAPA 2: Carregar PDBs existentes para comparacao
        # ==================================================================
        print("[05_rosettafold] Etapa 2: Carregando PDBs existentes...")

        existing_aso_atoms = None
        existing_sl_atoms = None

        if ASO_PDB.exists():
            try:
                existing_aso_atoms = parse_pdb(ASO_PDB)
                print(f"  MRL_ASO_001.pdb: {len(existing_aso_atoms)} atomos")
            except Exception as exc:
                warnings.append(f"Erro ao ler ASO PDB: {exc}")
                print(f"  AVISO: Nao foi possivel ler {ASO_PDB}: {exc}")
        else:
            print(f"  MRL_ASO_001.pdb nao encontrado (analise comparativa pulada)")

        if SL_RNA_PDB.exists():
            try:
                existing_sl_atoms = parse_pdb(SL_RNA_PDB)
                print(f"  sl_rna.pdb: {len(existing_sl_atoms)} atomos")
            except Exception as exc:
                warnings.append(f"Erro ao ler SL RNA PDB: {exc}")
                print(f"  AVISO: Nao foi possivel ler {SL_RNA_PDB}: {exc}")
        else:
            print(f"  sl_rna.pdb nao encontrado")

        # ==================================================================
        # ETAPA 3: Analise geometrica do duplex
        # ==================================================================
        print("[05_rosettafold] Etapa 3: Analisando geometria do duplex...")

        analysis = analyze_structure(
            duplex_atoms=duplex,
            free_rna_atoms=sl_rna_free,
            existing_pdb_atoms=existing_aso_atoms,
            config=config,
        )

        print(f"  Forma helicoidal: {analysis.helical_params.helix_form}")
        print(f"  Rise medio: {analysis.helical_params.mean_rise:.2f} A")
        print(f"  Twist medio: {analysis.helical_params.mean_twist:.1f} graus")
        print(f"  H-bonds (Watson-Crick): {analysis.n_watson_crick}")
        print(f"  H-bonds (total): {analysis.n_total_hbonds}")
        print(f"  RNase H acessivel: {analysis.rnase_h_accessible} "
              f"(score={analysis.rnase_h_score:.2f})")
        print(f"  Qualidade estrutural: {analysis.structural_quality}")

        # ==================================================================
        # ETAPA 4: Decomposicao de energia de ligacao
        # ==================================================================
        print("[05_rosettafold] Etapa 4: Decompondo energia de ligacao...")

        energy = decompose_binding_energy(config)

        print(f"  dG total: {energy.total_dg:.2f} kcal/mol")
        print(f"  dG experimental: {energy.experimental_dg:.2f} kcal/mol")
        print(f"  Desvio: {energy.dg_deviation:.1f}%")
        print(f"  Tm predita: {energy.predicted_tm:.1f} C")
        print(f"  Tm experimental: {energy.experimental_tm:.1f} C")
        print(f"  Classe de ligacao: {energy.binding_class}")
        print(f"  Funcional (dG < -15): {energy.is_functional}")

        # Componentes
        print("  Componentes de energia:")
        print(f"    Nearest-neighbor: {energy.nn_thermodynamic.value_kcal:.2f} kcal/mol")
        print(f"    Base stacking:    {energy.stacking_energy.value_kcal:.2f} kcal/mol")
        print(f"    H-bonds:          {energy.hbond_energy.value_kcal:.2f} kcal/mol")
        print(f"    Eletrostatica:    {energy.electrostatic_energy.value_kcal:.2f} kcal/mol")
        print(f"    Dessolvatacao:   +{energy.solvation_penalty.value_kcal:.2f} kcal/mol")
        print(f"    LNA estabiliz.:   {energy.lna_contribution.value_kcal:.2f} kcal/mol")
        print(f"    PS backbone:      {energy.ps_contribution.value_kcal:.2f} kcal/mol")

        # ==================================================================
        # ETAPA 5: Montar envelope de resultados
        # ==================================================================
        print("[05_rosettafold] Etapa 5: Montando resultados...")

        analysis_dict = _analysis_to_dict(analysis)
        energy_dict = _energy_to_dict(energy)

        envelope["status"] = "complete"
        envelope["artifacts"] = artifacts
        envelope["data"] = {
            "structural_analysis": analysis_dict,
            "binding_energy": energy_dict,
            "sequences": {
                "sl_rna": SL_RNA_SEQUENCE,
                "sl_rna_dna_template": SL_SEQUENCE,
                "aso": ASO_SEQUENCE,
                "target_region": SL_RNA_SEQUENCE[ASO_TARGET_START:ASO_TARGET_END],
                "target_start": ASO_TARGET_START,
                "target_end": ASO_TARGET_END,
            },
            "modifications": {
                "lna_positions": LNA_POSITIONS,
                "ps_backbone": True,
                "gapmer_design": {
                    "5prime_flank_lna": [0, 1, 2],
                    "dna_gap": list(range(3, 21)),
                    "3prime_flank_lna": [21, 22, 23],
                },
            },
            "existing_pdbs_analyzed": {
                "aso_pdb": str(ASO_PDB) if ASO_PDB.exists() else None,
                "sl_rna_pdb": str(SL_RNA_PDB) if SL_RNA_PDB.exists() else None,
            },
        }

        # Metricas de resumo
        envelope["metrics"] = {
            "n_atoms_duplex": len(duplex),
            "n_atoms_complex": len(complex_atoms),
            "n_hbonds_wc": analysis.n_watson_crick,
            "n_hbonds_total": analysis.n_total_hbonds,
            "helix_form": analysis.helical_params.helix_form,
            "dg_predicted": round(energy.total_dg, 2),
            "dg_experimental": energy.experimental_dg,
            "dg_deviation_pct": round(energy.dg_deviation, 1),
            "tm_predicted": round(energy.predicted_tm, 1),
            "tm_experimental": energy.experimental_tm,
            "rnase_h_accessible": analysis.rnase_h_accessible,
            "structural_quality": analysis.structural_quality,
            "binding_class": energy.binding_class,
        }

        # Conclusao
        envelope["summary"]["conclusion"] = _generate_conclusion(analysis, energy)
        envelope["summary"]["key_metrics"] = {
            "dg_predicted_kcal_mol": round(energy.total_dg, 2),
            "tm_predicted_celsius": round(energy.predicted_tm, 1),
            "helix_form": analysis.helical_params.helix_form,
            "rnase_h_accessible": analysis.rnase_h_accessible,
            "binding_class": energy.binding_class,
            "structural_quality": analysis.structural_quality,
        }

        envelope["warnings"] = warnings

    envelope["runtime_seconds"] = timer.elapsed

    # Salvar resultados JSON
    OUTPUT_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = write_result(envelope, output_dir=OUTPUT_RESULTS_DIR)
    print(f"\n[05_rosettafold] Resultado salvo em {output_path}")
    print(f"[05_rosettafold] Tempo total: {timer.elapsed:.2f}s")

    # Salvar tambem no diretorio padrao de resultados
    write_result(envelope)

    return envelope


def _generate_conclusion(
    analysis: StructureAnalysisResult,
    energy: BindingEnergyResult,
) -> str:
    """Gera conclusao textual a partir dos resultados da analise.

    Args:
        analysis: Resultado da analise estrutural.
        energy: Resultado da decomposicao de energia.

    Returns:
        Texto de conclusao para o envelope.
    """
    parts: list[str] = []

    # Geometria
    parts.append(
        f"O duplex ASO:SL RNA adota geometria {analysis.helical_params.helix_form}-form "
        f"com rise={analysis.helical_params.mean_rise:.2f} A e "
        f"twist={analysis.helical_params.mean_twist:.1f} graus, "
        f"consistente com hibrido RNA:DNA."
    )

    # H-bonds
    parts.append(
        f"Detectadas {analysis.n_watson_crick} pontes de hidrogenio Watson-Crick "
        f"e {analysis.n_total_hbonds} no total."
    )

    # Energia
    parts.append(
        f"A energia livre de ligacao predita (dG={energy.total_dg:.2f} kcal/mol) "
        f"e {energy.binding_class} e desvia {energy.dg_deviation:.1f}% do valor "
        f"experimental ({energy.experimental_dg} kcal/mol). "
        f"Tm predita: {energy.predicted_tm:.1f} C vs experimental {energy.experimental_tm} C."
    )

    # RNase H
    if analysis.rnase_h_accessible:
        parts.append(
            "O sulco menor do duplex e acessivel para a RNase H "
            f"(score={analysis.rnase_h_score:.2f}), indicando que o ASO "
            "pode mediar clivagem do SL RNA via mecanismo RNase H1."
        )
    else:
        parts.append(
            "AVISO: O sulco menor pode ter acessibilidade limitada para RNase H. "
            "Verificar experimentalmente."
        )

    # LNA
    n_lna = len(LNA_POSITIONS)
    parts.append(
        f"O design gapmer com {n_lna} posicoes LNA nos flancos estabiliza "
        f"o duplex em ~{abs(energy.lna_contribution.value_kcal):.1f} kcal/mol "
        "adicionais, enquanto o gap DNA central preserva a ativacao de RNase H."
    )

    return " ".join(parts)


if __name__ == "__main__":
    main()
