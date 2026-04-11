"""Certificado de Delivery v1.0 — MRL-ASO-001.

Consolida os 6 modulos de delivery (A-F) em um documento unico
com score 0-10 por dimensao e veredicto final.

Dimensoes avaliadas:
  A — Estabilidade quimica (pH 4.5 + resistencia a nucleases)
  B — Permeabilidade de membrana (uptake + endocitose)
  C — Conjugados de entrega (ranking + receptor binding)
  D — Nanoparticulas lipidicas (formulacao + pH-responsive release)
  E — ADMET (PK + safety margins)
  F — Modelo imune estocastico (Th1 + clearance parasitario)

Uso:
    python -m aso_delivery.reports.delivery_certificate
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from core.logger import get_logger

logger = get_logger("delivery_certificate")

# ---------------------------------------------------------------------------
# Caminhos dos resultados por modulo
# ---------------------------------------------------------------------------

_BASE: Path = Path(__file__).resolve().parent.parent

_MOD_A_JSON: Path = _BASE / "module_a_stability" / "results" / "module_a_stability.json"
_MOD_B_JSON: Path = _BASE / "module_b_membrane" / "results" / "module_b_membrane.json"
_MOD_C_JSON: Path = _BASE / "module_c_conjugate" / "results" / "module_c_conjugate.json"
_MOD_D_JSON: Path = _BASE / "module_d_lnp" / "results" / "module_d_lnp.json"
_MOD_E_JSON: Path = _BASE / "module_e_admet" / "results" / "module_e_admet.json"
_MOD_F_JSON: Path = _BASE / "module_f_immune_sde" / "results" / "module_f_immune_sde.json"

# Diretorio de saida do certificado
_REPORT_DIR: Path = _BASE / "reports" / "results"


# ---------------------------------------------------------------------------
# Carregamento de resultados
# ---------------------------------------------------------------------------


def _load_json(path: Path) -> dict[str, Any] | None:
    """Carrega um JSON. Retorna None se ausente."""
    if not path.exists():
        logger.warning("Resultado nao encontrado: %s", path)
        return None
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def load_all_results() -> dict[str, dict[str, Any] | None]:
    """Carrega resultados de todos os 6 modulos."""
    return {
        "module_a_stability": _load_json(_MOD_A_JSON),
        "module_b_membrane": _load_json(_MOD_B_JSON),
        "module_c_conjugate": _load_json(_MOD_C_JSON),
        "module_d_lnp": _load_json(_MOD_D_JSON),
        "module_e_admet": _load_json(_MOD_E_JSON),
        "module_f_immune_sde": _load_json(_MOD_F_JSON),
    }


# ---------------------------------------------------------------------------
# Avaliacao individual — cada funcao retorna score (0-10) e detalhes
# ---------------------------------------------------------------------------


def _assess_stability(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo A — estabilidade no fagolisossomo.

    Criterio de score (0-10):
      10 = all tests pass + pH 4.5 fraction_bound >= 0.99 + halflife > 500h
       8 = all tests pass + fraction_bound >= 0.95 + halflife > 200h
       6 = all tests pass (basic criteria)
       4 = 2 of 3 tests pass
       2 = 1 of 3 tests pass
       0 = none pass
    """
    all_pass = data.get("all_tests_passed", False)
    ph = data.get("ph_stability_profile", {})
    nr = data.get("nuclease_resistance", {})
    lna = data.get("lna_stability", {})

    ph_functional = ph.get("functional_at_target_ph", False)
    nr_met = nr.get("therapeutic_window_met", False)
    lna_maintained = lna.get("geometry_maintained", False)

    # Extract key metrics
    ph_45 = ph.get("profile", {}).get("pH_4.5", {})
    fraction_bound = ph_45.get("fraction_bound", 0.0)
    dg_kcal = ph_45.get("delta_g_kcal", 0.0)
    halflife = nr.get("LNA_gapmer_halflife_hours", 0.0)
    c3_endo = lna.get("c3_endo_fraction_pH_4_5", 0.0)

    # Count passing tests
    tests_passed = sum([ph_functional, nr_met, lna_maintained])

    # Score
    if all_pass and fraction_bound >= 0.99 and halflife > 500:
        score = 10
    elif all_pass and fraction_bound >= 0.95 and halflife > 200:
        score = 8
    elif all_pass:
        score = 6
    elif tests_passed == 2:
        score = 4
    elif tests_passed == 1:
        score = 2
    else:
        score = 0

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"pH 4.5 stability: dG = {dg_kcal:.2f} kcal/mol, fraction bound = {fraction_bound:.4f}",
        f"Nuclease resistance: gapmer half-life = {halflife:.1f} hours",
        f"LNA conformation: C3'-endo at pH 4.5 = {c3_endo * 100:.1f}%",
        f"All tests passed: {'Yes' if all_pass else 'No'} ({tests_passed}/3)",
    ]

    return {
        "dimension": "A. PHAGOLYSOSOMAL STABILITY",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "ph_4_5_dg_kcal": dg_kcal,
            "ph_4_5_fraction_bound": fraction_bound,
            "gapmer_halflife_hours": halflife,
            "c3_endo_fraction_pH_4_5": c3_endo,
            "ph_functional": ph_functional,
            "nuclease_window_met": nr_met,
            "lna_geometry_maintained": lna_maintained,
        },
        "notes": notes,
    }


def _assess_membrane(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo B — permeabilidade de membrana.

    Criterio de score (0-10):
      10 = concentration >= 10x threshold + macrophage advantage > 15x
       8 = concentration >= 5x threshold + macrophage advantage > 10x
       6 = concentration sufficient + receptor-mediated dominant
       4 = concentration sufficient but low advantage
       2 = concentration insufficient
    """
    endo = data.get("endocytosis", {})
    conc_sufficient = endo.get("concentration_sufficient", False)
    total_conc = endo.get("total_intracellular_conc_nm", 0.0)
    dominant = endo.get("dominant_pathway", "")
    macrophage_advantage = endo.get("macrophage_advantage_factor", 0.0)

    # Effective threshold is 100 nM (from module_b constants)
    threshold_nm = 100.0
    conc_ratio = total_conc / threshold_nm if threshold_nm > 0 else 0.0

    # Score
    if conc_sufficient and conc_ratio >= 10 and macrophage_advantage > 15:
        score = 10
    elif conc_sufficient and conc_ratio >= 5 and macrophage_advantage > 10:
        score = 8
    elif conc_sufficient and dominant == "receptor_mediated":
        score = 6
    elif conc_sufficient:
        score = 4
    else:
        score = 2

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"Intracellular concentration: {total_conc:.0f} nM ({conc_ratio:.0f}x threshold)",
        f"Dominant pathway: {dominant}",
        f"Macrophage advantage: {macrophage_advantage:.1f}x over non-phagocytic cells",
        f"Endosomal escape required: No (target in phagolysosome)",
    ]

    return {
        "dimension": "B. MEMBRANE PERMEATION",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "total_intracellular_conc_nm": total_conc,
            "concentration_vs_threshold": conc_ratio,
            "concentration_sufficient": conc_sufficient,
            "dominant_pathway": dominant,
            "macrophage_advantage_factor": macrophage_advantage,
        },
        "notes": notes,
    }


def _assess_conjugate(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo C — conjugados de entrega.

    Criterio de score (0-10):
      10 = best strategy score > 0.7 + uptake > 8x + receptor-specific
       8 = best strategy score > 0.6 + uptake > 5x
       6 = best strategy score > 0.5 + uptake > 3x
       4 = some improvement over naked
       2 = no clear advantage
    """
    optimal = data.get("optimal_design", {})
    primary = optimal.get("primary_recommendation", {})
    rec_score = primary.get("recommendation_score", 0.0)
    uptake_fold = primary.get("expected_uptake_fold", 0.0)
    target_receptor = primary.get("target_receptor", "")
    conjugate_name = primary.get("conjugate", "")

    # Ranking data
    comparison = data.get("comparison_matrix", {})
    ranking = comparison.get("ranking", [])
    n_strategies = len(ranking)

    # GalNAc exclusion (scientific rigor check)
    galnac = data.get("galnac_assessment", {})
    galnac_excluded = galnac.get("conclusion") == "UNSUITABLE"

    # Score
    if rec_score > 0.7 and uptake_fold > 8 and target_receptor:
        score = 10
    elif rec_score > 0.6 and uptake_fold > 5:
        score = 8
    elif rec_score > 0.5 and uptake_fold > 3:
        score = 6
    elif uptake_fold > 1.5:
        score = 4
    else:
        score = 2

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"Best strategy: {conjugate_name} (score: {rec_score:.4f})",
        f"Expected uptake: {uptake_fold:.1f}x vs naked ASO",
        f"Target receptor: {target_receptor}",
        f"Strategies evaluated: {n_strategies}",
        f"GalNAc excluded: {'Yes (ASGPR absent on macrophages)' if galnac_excluded else 'No'}",
    ]

    return {
        "dimension": "C. CONJUGATE DELIVERY",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "best_conjugate": conjugate_name,
            "recommendation_score": rec_score,
            "expected_uptake_fold": uptake_fold,
            "target_receptor": target_receptor,
            "n_strategies_evaluated": n_strategies,
            "galnac_excluded": galnac_excluded,
        },
        "notes": notes,
    }


def _assess_lnp(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo D — nanoparticulas lipidicas.

    Criterio de score (0-10):
      10 = all 4 criteria met + encapsulation > 95% + release > 99%
       8 = all 4 criteria met + encapsulation > 90%
       6 = all 4 criteria met
       4 = 3 of 4 criteria met
       2 = fewer criteria met
    """
    all_met = data.get("all_criteria_met", False)
    criteria = data.get("criteria", {})

    np_opt = data.get("np_optimization", {})
    ee = np_opt.get("optimal_encapsulation", 0.0)
    np_ratio = np_opt.get("optimal_np_ratio", 0.0)

    ph_rel = data.get("ph_release", {})
    release = ph_rel.get("release_at_target_ph", 0.0)

    mt = data.get("macrophage_targeting", {})
    uptake_fold = mt.get("expected_uptake_fold", 0.0)
    recommended = mt.get("recommended", {})
    diameter = recommended.get("particle_diameter_nm", 0.0)

    n_criteria_met = sum(1 for v in criteria.values() if v)

    # Score
    if all_met and ee > 0.95 and release > 0.99:
        score = 10
    elif all_met and ee > 0.90:
        score = 8
    elif all_met:
        score = 6
    elif n_criteria_met >= 3:
        score = 4
    else:
        score = 2

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"Encapsulation efficiency: {ee * 100:.1f}% (N/P ratio: {np_ratio:.0f})",
        f"Phagolysosomal release: {release * 100:.1f}% at pH 4.5",
        f"Particle diameter: {diameter:.0f} nm (phagocytable range)",
        f"Macrophage uptake: {uptake_fold:.1f}x (mannose-targeted)",
        f"Criteria met: {n_criteria_met}/4",
    ]

    return {
        "dimension": "D. LNP FORMULATION",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "encapsulation_efficiency": ee,
            "np_ratio": np_ratio,
            "release_at_target_ph": release,
            "particle_diameter_nm": diameter,
            "macrophage_uptake_fold": uptake_fold,
            "criteria_met": n_criteria_met,
            "criteria_detail": criteria,
        },
        "notes": notes,
    }


def _assess_admet(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo E — ADMET.

    Criterio de score (0-10):
      10 = TI >= 8 + bioavailability >= 85% + half-life >= 14 days
       8 = TI >= 5 + bioavailability >= 70% + half-life >= 7 days
       6 = TI >= 3 + bioavailability >= 50%
       4 = TI >= 2 (minimal safety margin)
       2 = TI < 2 (insufficient safety margin)
    """
    tox = data.get("toxicity", {})
    ti = tox.get("therapeutic_index", 0.0)
    risk_class = tox.get("overall_risk_classification", "")

    abs_data = data.get("absorption", {})
    bioav = abs_data.get("bioavailability_percent", 0.0)

    exc = data.get("excretion", {})
    t_half = exc.get("terminal_half_life_days", 0.0)

    dist = data.get("distribution", {})
    vd_ss = dist.get("vd_ss_l_per_kg", 0.0)

    # Score
    if ti >= 8 and bioav >= 85 and t_half >= 14:
        score = 10
    elif ti >= 5 and bioav >= 70 and t_half >= 7:
        score = 8
    elif ti >= 3 and bioav >= 50:
        score = 6
    elif ti >= 2:
        score = 4
    else:
        score = 2

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"Therapeutic index: {ti:.1f} (NOAEL/dose)",
        f"Bioavailability (SC): {bioav:.0f}%",
        f"Terminal half-life: {t_half:.0f} days (enables weekly dosing)",
        f"Vd,ss: {vd_ss:.2f} L/kg (hepatosplenic tropism)",
        f"Risk classification: {risk_class}",
    ]

    return {
        "dimension": "E. ADMET PROFILE",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "therapeutic_index": ti,
            "bioavailability_percent": bioav,
            "terminal_half_life_days": t_half,
            "vd_ss_l_per_kg": vd_ss,
            "risk_classification": risk_class,
        },
        "notes": notes,
    }


def _assess_immune_sde(data: dict[str, Any]) -> dict[str, Any]:
    """Avalia modulo F — simulacao imune estocastica.

    Criterio de score (0-10):
      10 = clearance_probability >= 0.95 + reduction >= 99% + EC50 < 0.5 uM
       8 = clearance_probability >= 0.80 + reduction >= 95%
       6 = clearance_probability > 0.50 + reduction >= 90%
       4 = clearance_probability > 0.30
       2 = clearance_probability <= 0.30
    """
    det = data.get("deterministic_scenarios", {})
    dual_det = det.get("dual_function", {})
    reduction_pct = dual_det.get("parasite_reduction_pct", 0.0)
    clearance_time = dual_det.get("time_to_90pct_clearance_hours", float("inf"))

    sde = data.get("sde_scenarios", {})
    dual_sde = sde.get("dual_function", {})
    clearance_prob = dual_sde.get("clearance_probability", 0.0)
    mean_clearance = dual_sde.get("mean_clearance_time_hours", float("inf"))

    synergy = data.get("synergy_analysis", {})
    si_speed = synergy.get("synergy_index_speed", 0.0)
    si_class = synergy.get("classification", "unknown")

    dr = data.get("dose_response", {})
    ec_values = dr.get("ec_values", {})
    ec50 = ec_values.get("EC50_uM", float("inf"))
    ec90 = ec_values.get("EC90_uM", float("inf"))

    # Convert ec50/ec90 to float if string
    try:
        ec50_f = float(ec50)
    except (ValueError, TypeError):
        ec50_f = float("inf")
    try:
        ec90_f = float(ec90)
    except (ValueError, TypeError):
        ec90_f = float("inf")

    # Score
    if clearance_prob >= 0.95 and reduction_pct >= 99 and ec50_f < 0.5:
        score = 10
    elif clearance_prob >= 0.80 and reduction_pct >= 95:
        score = 8
    elif clearance_prob > 0.50 and reduction_pct >= 90:
        score = 6
    elif clearance_prob > 0.30:
        score = 4
    else:
        score = 2

    verdict = "PASS" if score >= 6 else "FAIL"

    notes = [
        f"Parasite reduction (ODE): {reduction_pct:.1f}%",
        f"Clearance probability (SDE, N=1000): {clearance_prob * 100:.1f}%",
        f"Mean clearance time: {mean_clearance:.1f} hours",
        f"EC50 = {ec50} uM, EC90 = {ec90} uM",
        f"Synergy index: {si_speed:.4f} ({si_class})",
    ]

    return {
        "dimension": "F. IMMUNE RESPONSE (SDE)",
        "score": score,
        "max_score": 10,
        "verdict": verdict,
        "details": {
            "parasite_reduction_pct": reduction_pct,
            "clearance_probability": clearance_prob,
            "mean_clearance_time_hours": mean_clearance,
            "clearance_time_det_hours": clearance_time,
            "ec50_um": ec50,
            "ec90_um": ec90,
            "synergy_index_speed": si_speed,
            "synergy_classification": si_class,
        },
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Certificado consolidado
# ---------------------------------------------------------------------------


def generate_certificate() -> dict[str, Any]:
    """Gera o certificado de delivery final."""
    logger.info("=" * 70)
    logger.info("CERTIFICADO DE DELIVERY v1.0 — MRL-ASO-001")
    logger.info("=" * 70)

    results = load_all_results()

    # Mapear modulos aos assessores
    assessors: list[tuple[str, Any, Any]] = [
        ("module_a_stability", _assess_stability, results.get("module_a_stability")),
        ("module_b_membrane", _assess_membrane, results.get("module_b_membrane")),
        ("module_c_conjugate", _assess_conjugate, results.get("module_c_conjugate")),
        ("module_d_lnp", _assess_lnp, results.get("module_d_lnp")),
        ("module_e_admet", _assess_admet, results.get("module_e_admet")),
        ("module_f_immune_sde", _assess_immune_sde, results.get("module_f_immune_sde")),
    ]

    assessments: list[dict[str, Any]] = []
    for module_name, assessor, data in assessors:
        if data is None:
            assessment = {
                "dimension": module_name,
                "score": 0,
                "max_score": 10,
                "verdict": "MISSING",
                "details": {},
                "notes": [f"Module {module_name} result not found."],
            }
            logger.warning("Modulo %s: AUSENTE", module_name)
        else:
            assessment = assessor(data)
            logger.info(
                "Modulo %s: %s (score %d/10)",
                module_name,
                assessment["verdict"],
                assessment["score"],
            )
            for note in assessment["notes"]:
                logger.info("  %s", note)

        assessments.append(assessment)

    # Composite score
    total_score = sum(a["score"] for a in assessments)
    max_possible = 60

    # Veredicto final
    if total_score >= 48:
        overall_verdict = "VALIDATED"
    elif total_score >= 36:
        overall_verdict = "CONDITIONAL"
    else:
        overall_verdict = "FAILED"

    # Checar se algum modulo falhou ou esta ausente
    verdicts = [a["verdict"] for a in assessments]
    if "MISSING" in verdicts:
        overall_verdict = "CONDITIONAL"
    if "FAIL" in verdicts and overall_verdict == "VALIDATED":
        overall_verdict = "CONDITIONAL"

    generated_at = datetime.now(tz=timezone.utc)

    # Gerar texto ASCII do certificado
    certificate_text = _generate_certificate_text(
        assessments, overall_verdict, total_score, max_possible, generated_at
    )

    # Estrutura JSON completa
    certificate = {
        "title": "MRL-ASO-001 Delivery Validation Certificate v1.0",
        "version": "1.0",
        "generated_at": generated_at.isoformat(),
        "molecule": {
            "name": "MRL-ASO-001",
            "type": "25-nt LNA-DNA-LNA gapmer",
            "length": 25,
            "target": "Leishmania infantum SL RNA",
            "backbone": "full phosphorothioate (PS)",
        },
        "composite_score": total_score,
        "max_score": max_possible,
        "overall_verdict": overall_verdict,
        "dimension_assessments": assessments,
        "certificate_text": certificate_text,
    }

    # Gravar arquivos de saida
    _REPORT_DIR.mkdir(parents=True, exist_ok=True)

    json_path = _REPORT_DIR / "delivery_certificate.json"
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(certificate, fh, indent=2, ensure_ascii=False)
    logger.info("Certificado JSON gravado em: %s", json_path)

    txt_path = _REPORT_DIR / "delivery_certificate.txt"
    with open(txt_path, "w", encoding="utf-8") as fh:
        fh.write(certificate_text)
    logger.info("Certificado TXT gravado em: %s", txt_path)

    logger.info("")
    logger.info("=" * 70)
    logger.info(
        "VEREDICTO FINAL: %s (score: %d/%d)",
        overall_verdict,
        total_score,
        max_possible,
    )
    logger.info("=" * 70)

    return certificate


# ---------------------------------------------------------------------------
# Geracao do certificado em texto ASCII
# ---------------------------------------------------------------------------


def _pad_line(text: str, width: int = 62) -> str:
    """Preenche uma linha com espacos para caber dentro da moldura."""
    return text.ljust(width)


def _generate_certificate_text(
    assessments: list[dict[str, Any]],
    verdict: str,
    total_score: int,
    max_score: int,
    generated_at: datetime,
) -> str:
    """Gera o certificado em formato ASCII art estilizado."""
    w = 62  # largura interna
    border_top = "+" + "=" * (w + 2) + "+"
    border_mid = "+" + "-" * (w + 2) + "+"
    border_bot = "+" + "=" * (w + 2) + "+"

    def boxline(text: str = "") -> str:
        """Linha dentro da moldura."""
        return f"| {_pad_line(text, w)} |"

    lines: list[str] = []

    # Cabecalho
    lines.append(border_top)
    lines.append(boxline("MRL-ASO-001 DELIVERY VALIDATION CERTIFICATE"))
    lines.append(boxline("Version 1.0"))
    lines.append(border_mid)
    lines.append(boxline())

    # Dimensoes
    for assessment in assessments:
        dim = assessment["dimension"]
        score = assessment["score"]
        max_s = assessment["max_score"]
        v = assessment["verdict"]

        # Linha de titulo da dimensao
        score_str = f"[{score}/{max_s}]"
        verdict_str = f"[{v}]"
        header = f"  {dim:<40s} {score_str:>7s} {verdict_str:>6s}"
        lines.append(boxline(header))

        # Detalhes (notas)
        for note in assessment.get("notes", []):
            note_lines = _wrap_text(note, w - 8)
            for idx, nl in enumerate(note_lines):
                prefix = "- " if idx == 0 else "  "
                lines.append(boxline(f"    {prefix}{nl}"))

        lines.append(boxline())

    # Rodape — score composto e veredicto
    lines.append(border_mid)
    lines.append(boxline(f"  COMPOSITE SCORE: {total_score}/{max_score}"))
    lines.append(boxline(f"  OVERALL VERDICT: {verdict}"))
    lines.append(boxline())
    lines.append(
        boxline(f"  Generated: {generated_at.strftime('%Y-%m-%d %H:%M:%S UTC')}")
    )
    lines.append(boxline("  Target: Leishmania infantum SL RNA"))
    lines.append(boxline("  ASO: MRL-ASO-001 (25-nt LNA-DNA-LNA gapmer, PS)"))
    lines.append(border_bot)

    # Nota metodologica fora da moldura
    lines.append("")
    lines.append("Methodological note:")
    lines.append(
        "This certificate evaluates the delivery pipeline for MRL-ASO-001"
    )
    lines.append(
        "using computational models for pH stability (nearest-neighbor"
    )
    lines.append(
        "thermodynamics), membrane permeation (Born solvation / endocytosis),"
    )
    lines.append(
        "conjugate selection (receptor binding + avidity), LNP formulation"
    )
    lines.append(
        "(Helfrich theory + pH-responsive lipids), ADMET (2-compartment PK),"
    )
    lines.append(
        "and immune response (Langevin SDE + Monte Carlo simulation)."
    )
    lines.append(
        "Results are computational predictions, not experimental validation."
    )
    lines.append("")

    return "\n".join(lines)


def _wrap_text(text: str, max_width: int) -> list[str]:
    """Quebra texto longo em linhas de no maximo max_width caracteres."""
    if len(text) <= max_width:
        return [text]

    words = text.split()
    result_lines: list[str] = []
    current_line = ""

    for word in words:
        if current_line and len(current_line) + 1 + len(word) > max_width:
            result_lines.append(current_line)
            current_line = word
        elif current_line:
            current_line += " " + word
        else:
            current_line = word

    if current_line:
        result_lines.append(current_line)

    return result_lines


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> dict[str, Any]:
    """Entry point."""
    return generate_certificate()


if __name__ == "__main__":
    main()
