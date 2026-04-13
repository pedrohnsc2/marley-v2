"""Migracao unica: transforma os resultados estaticos atuais em 'Run Zero' (baseline).

Cria um run sintetico em results/runs/baseline/ com metadata.json
apontando para os dados existentes, e configura symlinks em results/latest/.

Uso:
    python scripts/migrate_to_runs.py
    python scripts/migrate_to_runs.py --dry-run

Este script e idempotente — pode ser executado varias vezes com seguranca.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
RUNS_DIR = RESULTS_DIR / "runs"
LATEST_DIR = RESULTS_DIR / "latest"
BASELINE_DIR = RUNS_DIR / "baseline"


def _create_baseline_metadata() -> dict:
    """Cria metadata sintetica para o run baseline (dados historicos)."""
    return {
        "run_id": "baseline",
        "pipeline": "all",
        "status": "completed",
        "created_at": "2026-04-01T00:00:00+00:00",
        "started_at": "2026-04-01T00:00:00+00:00",
        "completed_at": "2026-04-12T00:00:00+00:00",
        "git_sha": "historical",
        "parameters": {
            "organism": "leishmania_infantum",
            "species_name": "Leishmania infantum",
            "note": "Dados estaticos originais do pipeline — pre-parametrizacao",
        },
        "stages": [],
        "total_duration_s": 0.0,
        "output_dir": str(BASELINE_DIR / "outputs"),
        "notes": "Baseline: dados estaticos gerados durante o desenvolvimento do pipeline.",
        "tags": ["baseline", "historical", "leishmania_infantum"],
    }


def _setup_baseline_outputs(dry_run: bool = False) -> None:
    """Cria symlinks no baseline/outputs/ apontando para os dados reais."""
    outputs_dir = BASELINE_DIR / "outputs"

    if dry_run:
        print(f"  [DRY RUN] Would create: {outputs_dir}/")
        return

    outputs_dir.mkdir(parents=True, exist_ok=True)

    # Mapeia os diretorios de resultados existentes
    result_mappings = {
        # Dados no results/ raiz
        "construct": RESULTS_DIR / "construct",
        "docking": RESULTS_DIR / "docking",
        "aso": RESULTS_DIR / "aso",
        "rna": RESULTS_DIR / "rna",
        "filtered": RESULTS_DIR / "filtered",
        "conservation": RESULTS_DIR / "conservation",
        "immunogenicity": RESULTS_DIR / "immunogenicity",
        "report": RESULTS_DIR / "report",
    }

    # Cria symlinks para cada subdiretorio existente
    for name, source in result_mappings.items():
        link = outputs_dir / name
        if link.exists() or link.is_symlink():
            continue
        if source.exists():
            rel = os.path.relpath(source, outputs_dir)
            link.symlink_to(rel)
            print(f"  Linked: {name} -> {rel}")

    # CSVs soltos no results/ raiz
    for csv_file in RESULTS_DIR.glob("*.csv"):
        link = outputs_dir / csv_file.name
        if link.exists() or link.is_symlink():
            continue
        rel = os.path.relpath(csv_file, outputs_dir)
        link.symlink_to(rel)
        print(f"  Linked: {csv_file.name} -> {rel}")


def _setup_latest_symlinks(dry_run: bool = False) -> None:
    """Configura results/latest/ com symlinks para cada pipeline."""
    if dry_run:
        print(f"  [DRY RUN] Would create: {LATEST_DIR}/")
        return

    LATEST_DIR.mkdir(parents=True, exist_ok=True)

    # Para cada pipeline, aponta para os dados relevantes
    pipeline_mappings = {
        "vaccine": RESULTS_DIR,
        "drug": RESULTS_DIR,
        "docking": RESULTS_DIR,
        "rna": RESULTS_DIR,
        "aso_math": PROJECT_ROOT / "aso_math" / "results",
        "aso_delivery": PROJECT_ROOT / "aso_delivery" / "results",
        "platforms": PROJECT_ROOT / "vaccine_platforms" / "reports" / "results",
    }

    for name, target in pipeline_mappings.items():
        link = LATEST_DIR / name
        if link.is_symlink():
            link.unlink()
        elif link.exists():
            continue

        if target.exists():
            rel = os.path.relpath(target, LATEST_DIR)
            link.symlink_to(rel)
            print(f"  Latest: {name} -> {rel}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Migra resultados estaticos para runs/baseline/")
    parser.add_argument("--dry-run", action="store_true", help="Mostra o que seria feito")
    args = parser.parse_args()

    print("\n  Marley — Migracao para Run-Based Results")
    print("  " + "=" * 50)

    # 1. Baseline metadata
    print("\n  1. Criando baseline metadata...")
    if not args.dry_run:
        BASELINE_DIR.mkdir(parents=True, exist_ok=True)
        meta_path = BASELINE_DIR / "metadata.json"
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(_create_baseline_metadata(), fh, indent=2, ensure_ascii=False)
        print(f"     Escrito: {meta_path}")
    else:
        print(f"  [DRY RUN] Would write: {BASELINE_DIR / 'metadata.json'}")

    # 2. Baseline outputs (symlinks)
    print("\n  2. Linkando baseline outputs...")
    _setup_baseline_outputs(dry_run=args.dry_run)

    # 3. Latest symlinks
    print("\n  3. Configurando latest/ symlinks...")
    _setup_latest_symlinks(dry_run=args.dry_run)

    print("\n  Migracao concluida!")
    print("  " + "=" * 50)
    print()


if __name__ == "__main__":
    main()
