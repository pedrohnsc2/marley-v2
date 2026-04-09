I now have a thorough understanding of the entire codebase. Here is the complete architecture plan.

---

# Marley Platform -- Architecture Blueprint for Commercial ASO Validation

## 1. System Architecture

```
+===========================================================================+
|                          MARLEY PLATFORM                                   |
+===========================================================================+
|                                                                            |
|  CLIENTS                                                                   |
|  +-------------------+  +-------------------+  +---------------------+     |
|  | Web UI (Next.js)  |  | REST API Clients  |  | CLI (python -m)     |     |
|  | Upload sequence   |  | Programmatic use  |  | Local dev/research  |     |
|  | View reports      |  | POST /jobs        |  | python -m aso_math  |     |
|  | Download PDF      |  | GET /results      |  |   --target tcruzi   |     |
|  +--------+----------+  +--------+----------+  +----------+----------+     |
|           |                      |                         |               |
|  =========|======================|=========================|===========    |
|           v                      v                         v               |
|  +--------+----------------------+-------------------------+----------+    |
|  |                     API GATEWAY (FastAPI)                          |    |
|  |                                                                    |    |
|  |  /api/v1/jobs          - submit pipeline job                       |    |
|  |  /api/v1/jobs/{id}     - status + results                         |    |
|  |  /api/v1/organisms     - SL RNA registry CRUD                     |    |
|  |  /api/v1/benchmarks    - approved ASO benchmark results            |    |
|  |  /api/v1/certificates  - download PDF/JSON certificates            |    |
|  |  /api/v1/auth          - Supabase Auth proxy                       |    |
|  +---+------------------------+----------------------------+----------+    |
|      |                        |                            |               |
|      v                        v                            v               |
|  +---+----+           +------+------+           +---------+---------+      |
|  | JOB    |           | ASO_MATH    |           | REPORT            |      |
|  | QUEUE  |           | ENGINE      |           | GENERATOR         |      |
|  | (Redis |           | (Pure       |           | (WeasyPrint PDF)  |      |
|  |  or    |           |  Python,    |           |                   |      |
|  | Celery)|           |  stdlib)    |           | Jinja2 templates  |      |
|  +---+----+           +------+------+           | + Matplotlib viz  |      |
|      |                       |                   +---------+---------+      |
|      |                       |                             |               |
|      v                       v                             v               |
|  +---+----+           +------+------+           +---------+---------+      |
|  | WORKER |           | TARGET      |           | STORAGE           |      |
|  | POOL   +---------->+ REGISTRY    |           | (Supabase Storage |      |
|  |        |           | (SL RNA DB) |           |  + local /results)|      |
|  +--------+           +------+------+           +-------------------+      |
|                              |                                             |
|                              v                                             |
|                       +------+------+                                      |
|                       | SUPABASE    |                                      |
|                       | PostgreSQL  |                                      |
|                       |             |                                      |
|                       | - organisms |                                      |
|                       | - jobs      |                                      |
|                       | - results   |                                      |
|                       | - users     |                                      |
|                       | - benchmarks|                                      |
|                       +-------------+                                      |
+============================================================================+
```

### Data Flow for a Pipeline Run

```
User uploads target organism + sequence
        |
        v
API creates Job record (status=queued)
        |
        v
Worker picks up job, resolves TargetConfig:
   - SL RNA sequence from registry
   - Host transcriptome path
   - Organism-specific parameters (mutation rate, copy number, generation time)
        |
        v
aso_math engine runs 5 modules sequentially:
   01_thermodynamic_landscape(target_config) --> JSON
   02_selectivity_proof(target_config)       --> JSON
   03_evolutionary_conservation(target_config) --> JSON
   04_exhaustive_optimization(target_config)  --> JSON
   05_resistance_model(target_config)        --> JSON
        |
        v
Certificate generator reads 5 JSONs, produces:
   - math_certificate.json (machine-readable)
   - math_certificate.pdf  (branded, downloadable)
        |
        v
Results stored in Supabase + Storage
Job status updated to completed
User notified (webhook/email/polling)
```

---

## 2. Database Schema

All tables live in Supabase PostgreSQL. Row-Level Security (RLS) enforces user isolation for jobs/results. The organisms and benchmarks tables are public-read.

### Table: `organisms`
The central registry of SL RNA-bearing organisms.

```sql
CREATE TABLE organisms (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    slug            TEXT UNIQUE NOT NULL,       -- e.g. "leishmania_infantum"
    species_name    TEXT NOT NULL,              -- "Leishmania infantum"
    common_name     TEXT,                       -- "visceral leishmaniasis agent"
    taxonomic_group TEXT NOT NULL,              -- "kinetoplastida", "nematoda", "platyhelminthes"
    sl_rna_sequence TEXT NOT NULL,              -- 39 nt exon sequence
    sl_rna_length   INT NOT NULL,
    sl_rna_source   TEXT,                       -- "GenBank:X12345" or literature ref
    host_organism   TEXT,                       -- "Canis lupus familiaris", "Homo sapiens"
    disease_name    TEXT,                       -- "visceral leishmaniasis"
    affected_population_estimate BIGINT,        -- 2000000000 for all SL-bearing parasites
    mutation_rate   DOUBLE PRECISION NOT NULL DEFAULT 2.0e-9,  -- per nt per replication
    generation_time_hours DOUBLE PRECISION NOT NULL DEFAULT 12.0,
    sl_copy_number  INT NOT NULL DEFAULT 150,
    genome_gc_content DOUBLE PRECISION,
    transcriptome_available BOOLEAN DEFAULT FALSE,
    transcriptome_path TEXT,                    -- S3/Storage path if uploaded
    related_species TEXT[],                     -- array of organism slugs for conservation
    notes           TEXT,
    created_at      TIMESTAMPTZ DEFAULT now(),
    updated_at      TIMESTAMPTZ DEFAULT now()
);

-- Seed data examples:
-- L. infantum, L. donovani, L. major, L. braziliensis, L. mexicana, L. amazonensis
-- T. cruzi, T. brucei, T. vivax
-- C. fasciculata, L. seymouri
-- C. elegans (SL1 sequence: GGTTTAATTACCCAAGTTTGAG)
-- Ascaris suum, Brugia malayi (nematode SL sequences)
```

### Table: `benchmark_asos`
FDA-approved ASOs used to validate the thermodynamic model.

```sql
CREATE TABLE benchmark_asos (
    id                  UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    drug_name           TEXT NOT NULL UNIQUE,    -- "Nusinersen"
    brand_name          TEXT,                    -- "Spinraza"
    fda_approval_year   INT,
    target_gene         TEXT NOT NULL,           -- "SMN2 pre-mRNA"
    target_organism     TEXT DEFAULT 'Homo sapiens',
    aso_sequence        TEXT NOT NULL,
    aso_length          INT NOT NULL,
    chemistry           TEXT,                    -- "2'-MOE PS", "LNA gapmer", etc.
    mechanism           TEXT,                    -- "splice_switching", "rnase_h"
    published_tm        DOUBLE PRECISION,        -- from literature
    published_dg        DOUBLE PRECISION,
    predicted_tm        DOUBLE PRECISION,        -- our model output
    predicted_dg        DOUBLE PRECISION,
    tm_error            DOUBLE PRECISION,        -- |predicted - published|
    dg_error            DOUBLE PRECISION,
    benchmark_run_id    UUID,                    -- FK to jobs
    validation_status   TEXT DEFAULT 'pending',  -- pending, passed, failed
    references          JSONB,                   -- [{doi, title, year}]
    created_at          TIMESTAMPTZ DEFAULT now(),
    updated_at          TIMESTAMPTZ DEFAULT now()
);
```

### Table: `jobs`
Pipeline execution tracking.

```sql
CREATE TABLE jobs (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id         UUID REFERENCES auth.users(id),
    job_type        TEXT NOT NULL,              -- "full_validation", "benchmark", "single_module"
    status          TEXT NOT NULL DEFAULT 'queued',  -- queued, running, completed, failed
    organism_id     UUID REFERENCES organisms(id),
    organism_slug   TEXT,
    aso_sequence    TEXT,                       -- user-provided or auto-generated
    aso_name        TEXT,                       -- e.g. "MRL-ASO-001"
    target_start    INT,
    target_end      INT,
    config_override JSONB,                     -- override mutation_rate, etc.
    modules_requested TEXT[] DEFAULT ARRAY['01','02','03','04','05'],
    progress        JSONB DEFAULT '{}',        -- {"01": "completed", "02": "running", ...}
    results_path    TEXT,                       -- storage path for result artifacts
    certificate_path TEXT,                     -- path to PDF certificate
    score           DOUBLE PRECISION,          -- overall certificate score
    verdict         TEXT,                       -- VALIDATED, VALIDATED_WITH_NOTES, REQUIRES_REVIEW
    error_message   TEXT,
    started_at      TIMESTAMPTZ,
    completed_at    TIMESTAMPTZ,
    duration_seconds DOUBLE PRECISION,
    created_at      TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX idx_jobs_user_id ON jobs(user_id);
CREATE INDEX idx_jobs_status ON jobs(status);
```

### Table: `job_results`
Individual module outputs per job (one row per module per job).

```sql
CREATE TABLE job_results (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    job_id          UUID NOT NULL REFERENCES jobs(id) ON DELETE CASCADE,
    module_id       TEXT NOT NULL,              -- "01", "02", "03", "04", "05", "certificate"
    status          TEXT NOT NULL,              -- success, failed, skipped
    runtime_seconds DOUBLE PRECISION,
    result_json     JSONB NOT NULL,             -- full envelope from each module
    key_metrics     JSONB,                      -- extracted summary metrics
    verdict         TEXT,                       -- pass, flag, fail
    score           INT,
    created_at      TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX idx_job_results_job_id ON job_results(job_id);
```

### Table: `users` (managed by Supabase Auth, extended with profile)

```sql
CREATE TABLE user_profiles (
    id              UUID PRIMARY KEY REFERENCES auth.users(id),
    email           TEXT NOT NULL,
    full_name       TEXT,
    organization    TEXT,
    tier            TEXT DEFAULT 'free',        -- free, researcher, institution, enterprise
    monthly_jobs_used INT DEFAULT 0,
    monthly_jobs_limit INT DEFAULT 5,
    stripe_customer_id TEXT,
    stripe_subscription_id TEXT,
    created_at      TIMESTAMPTZ DEFAULT now()
);
```

### Table: `molecule_catalog`
Registry of validated ASO molecules across all organisms.

```sql
CREATE TABLE molecule_catalog (
    id              UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    molecule_name   TEXT NOT NULL,              -- "MRL-ASO-001"
    organism_id     UUID REFERENCES organisms(id),
    aso_sequence    TEXT NOT NULL,
    target_sequence TEXT NOT NULL,
    target_start    INT,
    target_end      INT,
    chemistry       TEXT,                       -- "LNA gapmer 4+17+4"
    dg_binding      DOUBLE PRECISION,
    tm_celsius      DOUBLE PRECISION,
    gc_content      DOUBLE PRECISION,
    certificate_score DOUBLE PRECISION,
    certificate_verdict TEXT,
    job_id          UUID REFERENCES jobs(id),
    is_public       BOOLEAN DEFAULT FALSE,
    created_at      TIMESTAMPTZ DEFAULT now()
);
```

---

## 3. API Design

Base URL: `/api/v1`

### Authentication
All endpoints except `GET /organisms` and `GET /benchmarks` require a Bearer token (Supabase JWT).

### Endpoints

```
POST   /auth/register              -- Supabase Auth signup
POST   /auth/login                 -- Supabase Auth signin
POST   /auth/refresh               -- Refresh JWT

GET    /organisms                  -- List all organisms in registry
GET    /organisms/{slug}           -- Get organism details + SL RNA
POST   /organisms                  -- Add organism (admin only)
PUT    /organisms/{slug}           -- Update organism (admin only)

POST   /jobs                       -- Submit a new pipeline job
GET    /jobs                       -- List user's jobs (paginated)
GET    /jobs/{id}                  -- Get job status + results
GET    /jobs/{id}/results          -- Get all module results for a job
GET    /jobs/{id}/results/{module} -- Get specific module result
DELETE /jobs/{id}                  -- Cancel queued job

GET    /jobs/{id}/certificate      -- Download PDF certificate
GET    /jobs/{id}/certificate.json -- Download JSON certificate

GET    /benchmarks                 -- List all benchmark ASOs and results
POST   /benchmarks/run             -- Trigger benchmark run (admin)
GET    /benchmarks/{drug_name}     -- Get specific benchmark result

GET    /molecules                  -- Browse public molecule catalog
GET    /molecules/{id}             -- Get molecule details

GET    /users/me                   -- Current user profile
PUT    /users/me                   -- Update profile
GET    /users/me/usage             -- Current billing period usage
```

### Job Submission Schema

```json
POST /jobs
{
    "organism_slug": "leishmania_infantum",
    "aso_sequence": "ACAGAAACTGATACTTATATAGCGT",      // optional: auto-design if omitted
    "aso_name": "MRL-ASO-001",                         // optional
    "target_start": 5,                                  // optional: auto-detect
    "target_end": 30,                                   // optional
    "modules": ["01", "02", "03", "04", "05"],         // optional: default all
    "config_override": {                                // optional
        "lna_5prime": 4,
        "lna_3prime": 4,
        "length_scan_min": 18,
        "length_scan_max": 27
    }
}
```

Response:
```json
{
    "job_id": "uuid",
    "status": "queued",
    "estimated_duration_seconds": 5,
    "poll_url": "/api/v1/jobs/{uuid}"
}
```

### Job Results Schema

```json
GET /jobs/{id}
{
    "job_id": "uuid",
    "status": "completed",
    "organism": "Leishmania infantum",
    "aso_name": "MRL-ASO-001",
    "score": 94.0,
    "verdict": "VALIDATED_WITH_NOTES",
    "duration_seconds": 4.48,
    "modules": {
        "01": {"status": "success", "verdict": "pass", "score": 90},
        "02": {"status": "success", "verdict": "pass", "score": 100},
        "03": {"status": "success", "verdict": "pass", "score": 95},
        "04": {"status": "success", "verdict": "pass", "score": 85},
        "05": {"status": "success", "verdict": "pass", "score": 100}
    },
    "certificate_url": "/api/v1/jobs/{uuid}/certificate",
    "created_at": "2026-04-09T14:00:00Z",
    "completed_at": "2026-04-09T14:00:04Z"
}
```

---

## 4. Generalization Strategy

The central problem: `aso_math/config.py` hardcodes L. infantum sequences, parameters, and file paths. The refactoring must preserve the stdlib-only, independently-runnable nature of each module while making them parameterizable.

### Step 1: Introduce `TargetConfig` dataclass

Create `/Users/pedronascimento/dev/marley/aso_math/target_config.py`:

```python
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

@dataclass(frozen=True)
class TargetConfig:
    """Immutable configuration for a single ASO validation run.

    Replaces all hardcoded values in config.py. A TargetConfig is
    constructed once at the start of a run and threaded through every
    module function.
    """
    # --- Organism identity ---
    organism_slug: str                     # "leishmania_infantum"
    species_name: str                      # "Leishmania infantum"

    # --- Target sequence ---
    sl_sequence: str                       # 39 nt SL RNA exon
    sl_length: int = 0                     # auto-computed from sl_sequence

    # --- ASO candidate ---
    aso_sequence: str = ""                 # to validate; empty = auto-design
    aso_name: str = ""
    aso_target_start: int = 0
    aso_target_end: int = 0

    # --- Organism-specific biology ---
    mutation_rate: float = 2.0e-9          # per nt per replication
    generation_time_hours: float = 12.0
    sl_copy_number: int = 150
    dg_functional_threshold: float = -15.0

    # --- Design parameters ---
    length_scan_min: int = 18
    length_scan_max: int = 27
    aso_concentration: float = 250e-9

    # --- LNA configuration ---
    lna_5prime: int = 4
    lna_3prime: int = 4

    # --- Cross-species conservation ---
    related_sl_sequences: dict[str, str] = field(default_factory=dict)
    divergence_time_mya: float = 350.0

    # --- File paths ---
    host_transcriptome_path: Path | None = None
    output_dir: Path = Path("results")

    # --- Known values for validation (optional) ---
    known_tm: float | None = None
    known_dg: float | None = None
    known_gc: float | None = None

    def __post_init__(self) -> None:
        if self.sl_length == 0:
            object.__setattr__(self, 'sl_length', len(self.sl_sequence))
        if self.aso_sequence and self.aso_target_end == 0:
            object.__setattr__(self, 'aso_target_end', self.sl_length)
```

### Step 2: Factory functions for known organisms

Create `/Users/pedronascimento/dev/marley/aso_math/organisms.py`:

```python
from aso_math.target_config import TargetConfig
from pathlib import Path

_REGISTRY: dict[str, dict] = {
    "leishmania_infantum": {
        "species_name": "Leishmania infantum",
        "sl_sequence": "AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG",
        "aso_sequence": "ACAGAAACTGATACTTATATAGCGT",
        "aso_name": "MRL-ASO-001",
        "aso_target_start": 5,
        "aso_target_end": 30,
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 12.0,
        "sl_copy_number": 150,
        "known_tm": 68.48,
        "known_dg": -27.97,
        "known_gc": 0.32,
        "related_sl_sequences": {
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
        },
    },
    "trypanosoma_cruzi": {
        "species_name": "Trypanosoma cruzi",
        "sl_sequence": "AACTAACGCTATTATTGATACAGTTTCTGTACTATATTG",
        "mutation_rate": 2.0e-9,
        "generation_time_hours": 24.0,
        "sl_copy_number": 200,
        "related_sl_sequences": { ... },  # same kinetoplastid set
    },
    "trypanosoma_brucei": { ... },
    "caenorhabditis_elegans": {
        "species_name": "Caenorhabditis elegans",
        "sl_sequence": "GGTTTAATTACCCAAGTTTGAG",  # SL1, 22 nt
        "mutation_rate": 1.7e-9,
        "generation_time_hours": 72.0,
        "sl_copy_number": 110,
        "related_sl_sequences": {
            "C. elegans SL1": "GGTTTAATTACCCAAGTTTGAG",
            "C. briggsae SL1": "GGTTTAATTACCCAAGTTTGAG",
            "Ascaris suum SL": "GGTTTAATTACCCAAGTTTGAG",
        },
    },
}

def get_target_config(slug: str, **overrides) -> TargetConfig:
    """Build a TargetConfig from the registry with optional overrides."""
    if slug not in _REGISTRY:
        raise ValueError(f"Unknown organism: {slug}. Available: {list(_REGISTRY.keys())}")
    params = {**_REGISTRY[slug], "organism_slug": slug}
    params.update(overrides)
    return TargetConfig(**params)

def list_organisms() -> list[str]:
    return list(_REGISTRY.keys())
```

### Step 3: Refactor each module's `main()` signature

Current pattern (every module):
```python
def main() -> dict[str, Any]:
    # reads from aso_math.config.*
```

New pattern:
```python
def main(config: TargetConfig | None = None) -> dict[str, Any]:
    if config is None:
        from aso_math.target_config import TargetConfig
        from aso_math.config import (SL_SEQUENCE, ASO_SEQUENCE, ...)
        config = TargetConfig(
            organism_slug="leishmania_infantum",
            species_name="Leishmania infantum",
            sl_sequence=SL_SEQUENCE,
            aso_sequence=ASO_SEQUENCE,
            ...
        )
    # Use config.sl_sequence instead of SL_SEQUENCE
    # Use config.aso_sequence instead of ASO_SEQUENCE
```

This preserves backward compatibility: `python -m aso_math.01_thermodynamic_landscape.run` still works with the L. infantum defaults. But programmatic callers and the API can inject a different `TargetConfig`.

### Step 4: Refactor `thermo.py` to accept parameters

The NN parameters (SantaLucia 1998) are universal constants -- they do not change per organism. What changes is `ASO_CONCENTRATION`. Refactor `compute_tm` and `compute_dg` to accept an optional `concentration` parameter:

```python
def compute_tm(sequence: str, concentration: float = 250e-9) -> float:
    ...

def compute_dg(sequence: str, temperature: float = 310.15) -> float:
    ...
```

### Step 5: Update `run_all.py` to accept `--organism` flag

```python
parser.add_argument("--organism", type=str, default="leishmania_infantum",
    help="Organism slug from registry (e.g., trypanosoma_cruzi)")
parser.add_argument("--sl-sequence", type=str, default=None,
    help="Override SL RNA sequence directly")
parser.add_argument("--aso-sequence", type=str, default=None,
    help="ASO sequence to validate (auto-designs if omitted)")
```

### Critical constraint preserved

`config.py` remains as-is -- it becomes the "L. infantum defaults" module. All NN_PARAMS, COMPLEMENT, R_GAS, etc. are universal constants that do not depend on organism. The organism-specific values (SL_SEQUENCE, MUTATION_RATE, etc.) get moved into `TargetConfig` while `config.py` provides the defaults for backward compatibility.

---

## 5. Benchmark Module Design

### Purpose
Validate the SantaLucia nearest-neighbor model against FDA-approved ASOs with published thermodynamic data. This is the single most critical credibility step -- if our model correctly predicts Tm/dG for Nusinersen, Mipomersen, Inotersen, and Volanesorsen, the predictions for MRL-ASO-001 carry real weight.

### File: `/Users/pedronascimento/dev/marley/aso_math/benchmark/run.py`

```python
"""Benchmark module: validate thermodynamic model against FDA-approved ASOs.

For each approved ASO with published Tm and dG values, we compute our
predicted values using the SantaLucia 1998 nearest-neighbor model and
report the error. The model is validated if mean absolute error is within
the known accuracy of NN predictions (~1-2 C for Tm, ~1-2 kcal/mol for dG).

Reference ASOs:
    Nusinersen (Spinraza) - 2'-MOE PS, splice-switching, SMN2
    Mipomersen (Kynamro) - 2'-MOE PS gapmer, RNase H, ApoB-100
    Inotersen (Tegsedi) - 2'-MOE PS gapmer, RNase H, TTR
    Volanesorsen (Waylivra) - 2'-MOE PS, RNase H, ApoC-III
"""

from dataclasses import dataclass

@dataclass
class BenchmarkASO:
    drug_name: str
    brand_name: str
    fda_year: int
    target_gene: str
    sequence: str               # DNA backbone sequence (without modifications)
    chemistry: str              # "2'-MOE PS", "LNA gapmer", etc.
    mechanism: str              # "splice_switching", "rnase_h"
    published_tm: float | None  # Celsius, from literature
    published_dg: float | None  # kcal/mol, from literature
    references: list[str]       # DOIs

APPROVED_ASOS: list[BenchmarkASO] = [
    BenchmarkASO(
        drug_name="Nusinersen",
        brand_name="Spinraza",
        fda_year=2016,
        target_gene="SMN2 pre-mRNA intron 7",
        sequence="TCACTTTCATAATGCTGG",  # 18-mer
        chemistry="2'-MOE phosphorothioate",
        mechanism="splice_switching",
        published_tm=65.0,   # Rigo et al. 2012
        published_dg=None,
        references=["10.1124/jpet.112.199679"],
    ),
    BenchmarkASO(
        drug_name="Mipomersen",
        brand_name="Kynamro",
        fda_year=2013,
        target_gene="ApoB-100 mRNA",
        sequence="GCCTCAGTCTGCTTCGCACC",  # 20-mer
        chemistry="2'-MOE PS gapmer (5-10-5)",
        mechanism="rnase_h",
        published_tm=None,
        published_dg=None,
        references=["10.1016/j.atherosclerosis.2007.06.026"],
    ),
    # Inotersen, Volanesorsen entries similarly structured
]

def run_benchmark() -> dict:
    """Run the thermodynamic model on all benchmark ASOs.

    Returns:
        Dictionary with per-ASO predictions and aggregate error metrics.
    """
    results = []
    for aso in APPROVED_ASOS:
        predicted_tm = compute_tm(aso.sequence)
        predicted_dg = compute_dg(aso.sequence)

        tm_error = abs(predicted_tm - aso.published_tm) if aso.published_tm else None
        dg_error = abs(predicted_dg - aso.published_dg) if aso.published_dg else None

        results.append({
            "drug_name": aso.drug_name,
            "sequence": aso.sequence,
            "length": len(aso.sequence),
            "predicted_tm": predicted_tm,
            "published_tm": aso.published_tm,
            "tm_error": tm_error,
            "predicted_dg": predicted_dg,
            "published_dg": aso.published_dg,
            "dg_error": dg_error,
        })

    # Aggregate metrics
    tm_errors = [r["tm_error"] for r in results if r["tm_error"] is not None]
    dg_errors = [r["dg_error"] for r in results if r["dg_error"] is not None]

    return {
        "benchmark_results": results,
        "mean_tm_error": sum(tm_errors) / len(tm_errors) if tm_errors else None,
        "mean_dg_error": sum(dg_errors) / len(dg_errors) if dg_errors else None,
        "model_validated": all(e < 3.0 for e in tm_errors) if tm_errors else False,
    }
```

### Important caveat for the benchmark

The SantaLucia model predicts DNA/DNA duplex thermodynamics. FDA-approved ASOs use chemical modifications (2'-MOE, PS backbone, LNA) that shift Tm and dG. The benchmark must account for this:

1. **Unmodified backbone comparison**: Run the NN model on the unmodified DNA sequence. Report the delta from published values. This delta is the "chemistry correction factor."
2. **Literature correction**: 2'-MOE raises Tm by approximately +1-2 C per modified nucleotide. PS backbone lowers Tm by approximately -0.5 C per linkage. LNA raises Tm by +3-8 C per substitution.
3. **Corrected prediction**: Apply the correction factor, report both raw and corrected predictions.

This becomes a table in the paper: "Table 1: Model validation against FDA-approved ASOs."

---

## 6. Multi-Organism Pipeline

### How it works concretely

The `TargetConfig` + `organisms.py` registry is all that is needed. The five modules already work on arbitrary sequences -- the only organism-specific inputs are:

| Parameter | L. infantum | T. cruzi | T. brucei | C. elegans |
|---|---|---|---|---|
| SL sequence | AACTAACGCTAT... (39 nt) | AACTAACGCTAT... (39 nt) | AACTAACGCTAT... (39 nt) | GGTTTAATTACC... (22 nt) |
| Mutation rate | 2.0e-9 | 2.0e-9 | 2.0e-9 | 1.7e-9 |
| Generation time | 12h | 24h | 6h | 72h |
| SL copy number | 150 | 200 | 200 | 110 |
| Related species | 11 kinetoplastids | 11 kinetoplastids | 11 kinetoplastids | 3 nematode SL1 |

### Module-specific adaptation

**Module 01 (Thermodynamic Landscape)**: No changes needed beyond accepting `config.sl_sequence` and `config.aso_sequence`. The mutation scan, double mutant scan, and length scan all work on any sequence.

**Module 02 (Selectivity)**: The `get_sl_rna_sequences()` function currently returns hardcoded sequences. Replace with `config.related_sl_sequences`. The human transcriptome screening works unchanged (the ASO is different, the human transcriptome is the same).

**Module 03 (Evolutionary Conservation)**: The `get_kinetoplastid_sl_sequences()` function must be generalized. For nematodes, it loads nematode SL sequences instead of kinetoplastid ones. The divergence time parameters come from `config.divergence_time_mya`.

**Module 04 (Exhaustive Optimization)**: Works on any sequence with `config.length_scan_min`, `config.length_scan_max`, and `config.sl_sequence`.

**Module 05 (Resistance Model)**: The most organism-specific module. Uses `config.mutation_rate`, `config.generation_time_hours`, `config.sl_copy_number`. The functional constraint logic (P(functional|mutation) = 0.0 for conserved positions) is driven by Module 03 results, not hardcoded.

### CLI invocation

```bash
# L. infantum (default, backward-compatible)
python -m aso_math.run_all

# T. cruzi
python -m aso_math.run_all --organism trypanosoma_cruzi

# Custom organism with inline sequence
python -m aso_math.run_all \
    --organism custom \
    --sl-sequence "GGTTTAATTACCCAAGTTTGAG" \
    --species-name "Caenorhabditis elegans"
```

---

## 7. PDF Report Generator

### Architecture

```
aso_math/reports/
    __init__.py
    math_certificate.py        -- existing JSON/TXT generator (keep)
    pdf_generator.py           -- new PDF generator
    templates/
        certificate.html       -- Jinja2 HTML template
        base.css               -- shared styles
    assets/
        marley_logo.svg
        watermark.svg
```

### Technology: WeasyPrint

WeasyPrint converts HTML+CSS to PDF. It handles page breaks, headers/footers, and @media print rules. It is a pure Python library that does not require headless Chrome.

### Template structure

The HTML template receives the same data dict that `math_certificate.py` already produces, plus rendered Matplotlib charts encoded as base64 PNG data URIs.

```python
# pdf_generator.py

from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import base64, io

def generate_pdf(certificate: dict, output_path: Path) -> Path:
    """Generate branded PDF certificate from validation results."""

    # 1. Render visualizations
    charts = {
        "heatmap": _render_heatmap(certificate),
        "length_scan": _render_length_scan(certificate),
        "entropy_profile": _render_entropy_profile(certificate),
        "pareto_front": _render_pareto_front(certificate),
        "resistance_timeline": _render_resistance_comparison(certificate),
    }

    # 2. Render HTML template
    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template("certificate.html")
    html_content = template.render(cert=certificate, charts=charts)

    # 3. Convert to PDF
    HTML(string=html_content).write_pdf(output_path)
    return output_path
```

### Visualizations (5 charts, one per module)

1. **Heatmap** (Module 01): 25x4 matrix of dG values, wildtype positions highlighted
2. **Length scan curve** (Module 01): dG vs. length, MRL-ASO-001 position marked
3. **Entropy profile** (Module 03): Shannon entropy per position, target region shaded
4. **Pareto front** (Module 04): scatter plot of dG vs Tm, Pareto-optimal points connected
5. **Resistance timeline** (Module 05): bar chart comparing time-to-resistance across drugs

---

## 8. Complete File Structure

```
marley/
|-- aso_math/                          # CORE ENGINE (stdlib-only, ~5000 LOC)
|   |-- __init__.py
|   |-- config.py                      # Universal constants (NN params, R, T, bases)
|   |-- target_config.py              # NEW: TargetConfig dataclass
|   |-- organisms.py                   # NEW: organism registry + factory
|   |-- thermo.py                      # Thermodynamic calculations (parameterized)
|   |-- envelope.py                    # JSON envelope utilities
|   |-- run_all.py                     # Orchestrator (accepts --organism)
|   |-- 01_thermodynamic_landscape/
|   |   |-- __init__.py
|   |   |-- run.py                     # main(config: TargetConfig | None)
|   |-- 02_selectivity_proof/
|   |   |-- __init__.py
|   |   |-- run.py
|   |-- 03_evolutionary_conservation/
|   |   |-- __init__.py
|   |   |-- run.py
|   |-- 04_exhaustive_optimization/
|   |   |-- __init__.py
|   |   |-- run.py
|   |-- 05_resistance_model/
|   |   |-- __init__.py
|   |   |-- run.py
|   |-- benchmark/                     # NEW: FDA-approved ASO benchmarks
|   |   |-- __init__.py
|   |   |-- approved_asos.py           # Nusinersen, Mipomersen, etc.
|   |   |-- run.py                     # run_benchmark() -> dict
|   |   |-- corrections.py            # Chemistry correction factors
|   |-- reports/
|   |   |-- __init__.py
|   |   |-- math_certificate.py        # Existing JSON/TXT certificate
|   |   |-- pdf_generator.py          # NEW: WeasyPrint PDF
|   |   |-- templates/
|   |   |   |-- certificate.html
|   |   |   |-- base.css
|   |   |-- assets/
|   |       |-- marley_logo.svg
|   |-- results/                       # Module outputs (gitignored)
|
|-- api/                               # NEW: FastAPI backend
|   |-- __init__.py
|   |-- main.py                        # FastAPI app, CORS, middleware
|   |-- config.py                      # API settings (env vars)
|   |-- auth.py                        # Supabase JWT validation
|   |-- routers/
|   |   |-- __init__.py
|   |   |-- jobs.py                    # POST/GET /jobs endpoints
|   |   |-- organisms.py              # GET /organisms endpoints
|   |   |-- benchmarks.py            # GET /benchmarks endpoints
|   |   |-- certificates.py          # GET /certificates (PDF download)
|   |   |-- molecules.py             # GET /molecules catalog
|   |   |-- users.py                  # GET/PUT /users/me
|   |-- workers/
|   |   |-- __init__.py
|   |   |-- pipeline_worker.py        # Celery/ARQ task that runs aso_math
|   |   |-- celery_app.py            # Celery configuration
|   |-- schemas/
|   |   |-- __init__.py
|   |   |-- job.py                    # Pydantic models for request/response
|   |   |-- organism.py
|   |   |-- benchmark.py
|   |   |-- user.py
|   |-- db/
|       |-- __init__.py
|       |-- supabase_client.py        # Async Supabase client
|       |-- migrations/               # SQL migration files
|           |-- 001_organisms.sql
|           |-- 002_benchmark_asos.sql
|           |-- 003_jobs.sql
|           |-- 004_user_profiles.sql
|           |-- 005_molecule_catalog.sql
|
|-- web/                               # Next.js frontend (EXISTING, expanded)
|   |-- app/
|   |   |-- page.tsx                   # Landing/dashboard
|   |   |-- (auth)/
|   |   |   |-- login/page.tsx
|   |   |   |-- register/page.tsx
|   |   |-- dashboard/
|   |   |   |-- page.tsx              # User dashboard
|   |   |   |-- jobs/
|   |   |   |   |-- page.tsx          # Job list
|   |   |   |   |-- [id]/page.tsx     # Job detail + results
|   |   |   |-- new/page.tsx          # New job wizard
|   |   |-- organisms/
|   |   |   |-- page.tsx              # Browse organisms
|   |   |-- benchmarks/
|   |   |   |-- page.tsx              # Benchmark results
|   |   |-- pricing/
|   |   |   |-- page.tsx
|   |   |-- api/
|   |       |-- data/route.ts         # Existing local file reader
|   |-- components/
|   |   |-- nav.tsx
|   |   |-- kpi-card.tsx
|   |   |-- sequence-input.tsx        # NEW: FASTA/raw sequence input
|   |   |-- job-status-badge.tsx      # NEW
|   |   |-- heatmap-chart.tsx         # NEW: interactive heatmap
|   |   |-- entropy-chart.tsx         # NEW: entropy profile
|   |   |-- pareto-chart.tsx          # NEW: Pareto front scatter
|   |-- lib/
|       |-- supabase.ts
|       |-- types.ts
|       |-- api-client.ts            # NEW: typed fetch wrapper for /api/v1
|
|-- core/                              # EXISTING shared utilities
|   |-- __init__.py
|   |-- logger.py
|   |-- models.py
|   |-- db.py
|   |-- uniprot.py
|   |-- codon_tables.py
|   |-- structure.py
|
|-- rna_entropy/                       # EXISTING ASO design pipeline
|-- pipeline/                          # EXISTING vaccine pipeline
|-- drug_targets/                      # EXISTING drug target pipeline
|-- orchestrator/                      # EXISTING CI/CD scripts
|
|-- tests/                             # EXPANDED test suite
|   |-- __init__.py
|   |-- conftest.py                   # Shared fixtures
|   |-- aso_math/
|   |   |-- __init__.py
|   |   |-- test_thermo.py            # NN model unit tests
|   |   |-- test_config.py
|   |   |-- test_target_config.py
|   |   |-- test_organisms.py
|   |   |-- test_01_thermodynamic.py
|   |   |-- test_02_selectivity.py
|   |   |-- test_03_conservation.py
|   |   |-- test_04_optimization.py
|   |   |-- test_05_resistance.py
|   |   |-- test_benchmark.py
|   |   |-- test_certificate.py
|   |   |-- test_pdf_generator.py
|   |-- api/
|   |   |-- test_jobs.py
|   |   |-- test_organisms.py
|   |   |-- test_auth.py
|   |-- integration/
|       |-- test_full_pipeline.py     # End-to-end: submit job, get PDF
|       |-- test_multi_organism.py
|
|-- data/                              # Data files (gitignored for large ones)
|   |-- organisms/                    # NEW: per-organism data
|   |   |-- leishmania_infantum/
|   |   |   |-- transcriptome.fasta
|   |   |-- trypanosoma_cruzi/
|   |   |-- trypanosoma_brucei/
|   |   |-- caenorhabditis_elegans/
|   |-- benchmarks/                   # NEW: published ASO data
|       |-- nusinersen.json
|       |-- mipomersen.json
|
|-- docs/                              # EXISTING documentation
|-- requirements.txt                   # EXISTING + new deps
|-- requirements-api.txt              # NEW: API-specific deps
|-- pyproject.toml                    # NEW: replace requirements.txt
|-- Dockerfile                        # NEW: API + worker container
|-- docker-compose.yml                # NEW: full stack local dev
|-- .env.example
|-- README.md
```

---

## 9. Technology Choices with Justification

| Layer | Choice | Justification |
|---|---|---|
| **Core engine** | Python 3.14, stdlib only | Already proven at 4.48s. No heavy deps means zero supply-chain risk for the scientific core. NN params are math, not ML. |
| **API framework** | FastAPI | Async, auto-generated OpenAPI docs, Pydantic validation, mainstream in biotech SaaS. FastAPI handles file uploads (FASTA) and streaming responses (SSE for job progress) natively. |
| **Task queue** | Celery + Redis | Pipeline runs take ~5s but must not block HTTP. Celery is battle-tested. Redis is already common in Supabase stacks. For initial simplicity, ARQ (async Redis queue) is a lighter alternative. |
| **Database** | Supabase PostgreSQL | Already integrated. RLS gives per-user data isolation for free. PostgREST API means the web frontend can query directly for read-heavy pages. |
| **Auth** | Supabase Auth | Already configured. Supports email/password, OAuth, magic links. JWT validation in FastAPI via `supabase.auth.get_user()`. |
| **PDF generation** | WeasyPrint | Pure Python, CSS-based layout, handles pages/headers/footers. No headless browser needed. Alternative: reportlab (more control, more code). |
| **Visualizations** | Matplotlib (PDF) + Recharts (web) | Matplotlib for server-side static charts in PDFs. Recharts for interactive web charts. Both are standard in their ecosystems. |
| **Frontend** | Next.js 14 + Tailwind | Already exists. App Router for SSR/RSC, Tailwind for rapid styling. No framework migration needed. |
| **Payments** | Stripe | Industry standard. Checkout Sessions for one-time purchases, Subscriptions for tiers. Webhooks update `user_profiles.tier`. |
| **Hosting** | Vercel (web) + Railway/Fly.io (API) | Vercel for Next.js is zero-config. Railway or Fly.io for the FastAPI + Celery workers with Docker. Supabase is already cloud-hosted. |
| **Testing** | pytest + pytest-cov + httpx | pytest is standard. httpx for async API testing. Coverage target >80% for aso_math. |

---

## 10. Implementation Phases

### Phase 1: Scientific Credibility (Weeks 1-3)
**Goal**: Produce benchmark data and multi-organism runs before any publication or patent filing.

| Week | Task | Deliverable | Est. Effort |
|---|---|---|---|
| 1 | Create `aso_math/benchmark/approved_asos.py` with Nusinersen, Mipomersen, Inotersen, Volanesorsen sequences and published Tm/dG | Benchmark data file | 4h |
| 1 | Create `aso_math/benchmark/run.py` that runs `compute_tm`/`compute_dg` on approved ASOs and computes error | Benchmark runner + JSON output | 6h |
| 1 | Research published Tm/dG values from literature for each approved ASO; document DOIs | Annotated reference table | 8h |
| 1 | Implement chemistry correction factors (2'-MOE, PS, LNA adjustments) | `corrections.py` | 4h |
| 2 | Create `TargetConfig` dataclass and `organisms.py` registry | `target_config.py`, `organisms.py` | 4h |
| 2 | Refactor Module 01 to accept `TargetConfig` (backward-compatible) | Updated `01_*/run.py` | 3h |
| 2 | Refactor Modules 02-05 to accept `TargetConfig` | Updated run.py for each | 8h |
| 2 | Run full pipeline for T. cruzi, validate results make biological sense | T. cruzi results JSON | 4h |
| 3 | Run full pipeline for T. brucei | T. brucei results JSON | 2h |
| 3 | Run pipeline for C. elegans SL1 (nematode demonstration) | C. elegans results JSON | 4h |
| 3 | Compile benchmark + multi-organism results into patent-ready format | Summary document | 8h |

**Dependencies**: None -- this phase uses only existing code + refactoring.
**Deliverables**: Benchmark validation table, 4 organism pipeline runs, patent claims structure.

### Phase 2: Software Generalization + Testing (Weeks 4-6)
**Goal**: Production-quality codebase with test coverage.

| Week | Task | Deliverable | Est. Effort |
|---|---|---|---|
| 4 | Write unit tests for `thermo.py` (compute_tm, compute_dg, edge cases) | `tests/aso_math/test_thermo.py` | 6h |
| 4 | Write unit tests for each module (known-answer tests) | 5 test files | 12h |
| 4 | Write tests for `TargetConfig` and `organisms.py` | 2 test files | 3h |
| 5 | Write integration test: full pipeline L. infantum -> certificate | `test_full_pipeline.py` | 4h |
| 5 | Write integration test: multi-organism runs | `test_multi_organism.py` | 4h |
| 5 | Set up CI (GitHub Actions): lint + test + coverage report | `.github/workflows/ci.yml` | 4h |
| 5 | Expand organism registry to 15+ species with literature-sourced SL sequences | Updated `organisms.py` | 8h |
| 6 | Add `pyproject.toml`, consolidate dependency management | `pyproject.toml` | 2h |
| 6 | Refactor `config.py`: separate universal constants from organism defaults cleanly | Updated `config.py` | 3h |
| 6 | Code review, docstring audit, type annotation completeness | Clean codebase | 6h |

**Dependencies**: Phase 1 complete (TargetConfig exists).
**Deliverables**: >80% test coverage, CI green, 15+ organisms in registry.

### Phase 3: Product Interface (Weeks 7-10)
**Goal**: Working web application where a user can submit a sequence and get a PDF.

| Week | Task | Deliverable | Est. Effort |
|---|---|---|---|
| 7 | FastAPI skeleton: project structure, config, Supabase client | `api/` directory | 6h |
| 7 | Implement `POST /jobs` and `GET /jobs/{id}` endpoints | Job CRUD | 8h |
| 7 | Implement Celery worker that runs `aso_math.run_all` | `pipeline_worker.py` | 6h |
| 8 | Implement `GET /organisms` endpoint backed by Supabase | Organisms API | 3h |
| 8 | Implement `GET /benchmarks` endpoint | Benchmarks API | 3h |
| 8 | Build PDF generator with WeasyPrint + Jinja2 template | `pdf_generator.py` + template | 12h |
| 8 | Implement `GET /jobs/{id}/certificate` (PDF download) | Certificate endpoint | 3h |
| 9 | Web: new job wizard page (organism selector, sequence input, submit) | `web/app/dashboard/new/page.tsx` | 10h |
| 9 | Web: job detail page (status, module results, download certificate) | `web/app/dashboard/jobs/[id]/page.tsx` | 10h |
| 9 | Web: interactive heatmap chart component (Recharts) | `heatmap-chart.tsx` | 6h |
| 10 | Web: entropy profile + Pareto front charts | 2 chart components | 8h |
| 10 | API tests: job submission, status polling, PDF download | `tests/api/` | 6h |
| 10 | Docker Compose for local full-stack development | `docker-compose.yml` | 4h |

**Dependencies**: Phase 2 complete (tests, generalized engine).
**Deliverables**: Working API + web UI, PDF certificate download, Docker setup.

### Phase 4: Commercial Infrastructure (Weeks 11-14)
**Goal**: Monetizable product.

| Week | Task | Deliverable | Est. Effort |
|---|---|---|---|
| 11 | Supabase Auth integration (register, login, JWT in API) | Auth flow | 8h |
| 11 | User profile + usage metering (track jobs/month) | `user_profiles` table + API | 6h |
| 12 | Stripe integration: pricing tiers, Checkout Sessions | Stripe webhooks | 12h |
| 12 | Rate limiting + tier enforcement in API | Middleware | 4h |
| 13 | Landing page with pricing table | `web/app/pricing/page.tsx` | 8h |
| 13 | Organisms browse page (public) | `web/app/organisms/page.tsx` | 4h |
| 13 | Benchmarks results page (public, for credibility) | `web/app/benchmarks/page.tsx` | 4h |
| 14 | Production deployment: Vercel + Railway + domain + SSL | Live site | 8h |
| 14 | Monitoring: error tracking (Sentry), uptime (Better Uptime) | Observability | 4h |

**Dependencies**: Phase 3 complete (working web app).
**Deliverables**: Live product at marley.bio (or similar), payment processing, user accounts.

### Phase 5: Scientific Publication (Weeks 15-20)
**Goal**: Peer-reviewed paper.

| Week | Task | Deliverable | Est. Effort |
|---|---|---|---|
| 15-16 | Write manuscript: Methods + Results (benchmark table, multi-organism, MRL-ASO-001) | Draft sections | 20h |
| 17 | Figures: publication-quality versions of the 5 module visualizations | 5 figures | 12h |
| 18 | Write Introduction + Discussion | Draft sections | 12h |
| 19 | Internal review, address gaps, run additional analyses if needed | Revised draft | 10h |
| 20 | Submit to PLoS NTD or Nucleic Acids Research | Submitted manuscript | 4h |

**Dependencies**: Phase 1 complete (benchmark + multi-organism data).

---

## 11. Risk Assessment

### R1: SantaLucia model inaccuracy for modified nucleotides
**Severity**: HIGH
**Likelihood**: MEDIUM
**Description**: The NN model is for DNA/DNA duplexes. ASOs use 2'-MOE, PS, LNA modifications that change thermodynamics. If benchmark errors exceed 5 C for Tm, the model loses credibility.
**Mitigation**: (a) Report unmodified predictions alongside correction-adjusted predictions. (b) Use published correction factors from McTigue 2004 (LNA) and Freier & Altmann 1997 (2'-MOE). (c) Frame the benchmark honestly: "Our model predicts unmodified backbone thermodynamics within X C; correction factors for Y chemistry bring predictions within Z C of published values." (d) If errors are large, implement Allawi & SantaLucia (1997) DNA/RNA hybrid parameters instead of DNA/DNA.

### R2: Lack of published Tm/dG for all approved ASOs
**Severity**: MEDIUM
**Likelihood**: HIGH
**Description**: Not all FDA-approved ASOs have publicly available Tm and dG measurements. Nusinersen has good literature data; others may require digging or may be proprietary.
**Mitigation**: (a) Prioritize Nusinersen as the primary validation target (best literature coverage). (b) Use Mipomersen and Inotersen as secondary (Ionis/ISIS has published some values). (c) If published values are unavailable, compute predicted values and compare against known clinical efficacy instead. (d) Document clearly which values are literature-sourced vs. predicted.

### R3: Nematode SL RNA is structurally different from kinetoplastid SL RNA
**Severity**: LOW
**Likelihood**: CERTAIN
**Description**: The C. elegans SL1 sequence (22 nt) is shorter and unrelated to the trypanosomatid SL (39 nt). Module 03 (conservation analysis) needs different reference species. Module 05 (resistance model) needs nematode-specific parameters.
**Mitigation**: The `TargetConfig` + `organisms.py` architecture handles this cleanly. Each organism brings its own related species set and biology parameters. Module 03 computes conservation among whatever species are provided; it does not assume kinetoplastid phylogeny.

### R4: Performance at scale (many concurrent jobs)
**Severity**: LOW
**Likelihood**: LOW (initially)
**Description**: Each pipeline run takes ~5s on a single core. At 100 concurrent users, that is 100 cores-worth of compute.
**Mitigation**: (a) Celery workers scale horizontally. (b) Each job is independent and CPU-bound (no I/O bottleneck). (c) At $X/month per tier, even modest revenue pays for compute. (d) Add caching: identical (organism, sequence) pairs return cached results instantly.

### R5: IP risk -- someone clones the pipeline
**Severity**: MEDIUM
**Likelihood**: MEDIUM
**Description**: The core algorithm (SantaLucia NN model) is published science. The value is in the 5-module validation framework, organism registry, and benchmark data -- but these could be replicated.
**Mitigation**: (a) File provisional patent before publication (Phase 1 deliverable). Claims cover the specific 5-proof validation methodology and the multi-barrier resistance model. (b) Speed to market: first-mover advantage in SL RNA ASO design tools. (c) Network effects: the organism registry and benchmark database grow with users. (d) The paper itself becomes a citation barrier.

### R6: Supabase vendor lock-in
**Severity**: LOW
**Likelihood**: LOW
**Description**: Core DB, auth, and storage depend on Supabase.
**Mitigation**: (a) All SQL is standard PostgreSQL -- portable to any Postgres host. (b) Auth can be replaced with any JWT provider. (c) Storage is S3-compatible. (d) The aso_math engine has zero Supabase dependency (stdlib only).

### R7: WeasyPrint rendering quality
**Severity**: LOW
**Likelihood**: MEDIUM
**Description**: Complex CSS (multi-column layouts, SVG charts) can render differently across WeasyPrint versions.
**Mitigation**: (a) Keep the certificate template simple: single-column, sequential sections. (b) Charts are pre-rendered as PNG (Matplotlib), not inline SVG. (c) Pin WeasyPrint version. (d) Fallback: use reportlab for pixel-perfect control if WeasyPrint proves unreliable.

### R8: Regulatory and clinical liability
**Severity**: HIGH
**Likelihood**: LOW
**Description**: Users might interpret certificate scores as clinical endorsements.
**Mitigation**: (a) Every certificate includes a mandatory disclaimer: "Computational predictions only. Not validated experimentally. Not medical advice." (b) Terms of Service disclaim clinical use. (c) The scoring is explicitly documented as a relative ranking tool, not an absolute efficacy predictor.

---

## Summary of Key Files to Create/Modify

**New files** (in priority order):
1. `/Users/pedronascimento/dev/marley/aso_math/target_config.py` -- TargetConfig dataclass
2. `/Users/pedronascimento/dev/marley/aso_math/organisms.py` -- organism registry
3. `/Users/pedronascimento/dev/marley/aso_math/benchmark/run.py` -- benchmark against approved ASOs
4. `/Users/pedronascimento/dev/marley/aso_math/benchmark/approved_asos.py` -- benchmark data
5. `/Users/pedronascimento/dev/marley/aso_math/benchmark/corrections.py` -- chemistry corrections
6. `/Users/pedronascimento/dev/marley/aso_math/reports/pdf_generator.py` -- WeasyPrint PDF
7. `/Users/pedronascimento/dev/marley/api/main.py` -- FastAPI application
8. `/Users/pedronascimento/dev/marley/api/routers/jobs.py` -- Job endpoints
9. `/Users/pedronascimento/dev/marley/api/workers/pipeline_worker.py` -- Celery worker

**Files to modify** (generalization):
1. `/Users/pedronascimento/dev/marley/aso_math/01_thermodynamic_landscape/run.py` -- accept TargetConfig
2. `/Users/pedronascimento/dev/marley/aso_math/02_selectivity_proof/run.py` -- accept TargetConfig
3. `/Users/pedronascimento/dev/marley/aso_math/03_evolutionary_conservation/run.py` -- accept TargetConfig
4. `/Users/pedronascimento/dev/marley/aso_math/04_exhaustive_optimization/run.py` -- accept TargetConfig
5. `/Users/pedronascimento/dev/marley/aso_math/05_resistance_model/run.py` -- accept TargetConfig
6. `/Users/pedronascimento/dev/marley/aso_math/run_all.py` -- --organism flag
7. `/Users/pedronascimento/dev/marley/aso_math/thermo.py` -- parameterize concentration
8. `/Users/pedronascimento/dev/marley/aso_math/reports/math_certificate.py` -- accept TargetConfig

**Files that remain unchanged**:
- `/Users/pedronascimento/dev/marley/aso_math/config.py` -- becomes "defaults" module, constants are universal
- `/Users/pedronascimento/dev/marley/aso_math/envelope.py` -- already generic
- `/Users/pedronascimento/dev/marley/core/logger.py` -- already generic

The generalization strategy is designed so that zero existing functionality breaks. Every `main()` function defaults to L. infantum when called without arguments. The `TargetConfig` is an additive parameter, not a breaking change. This means the scientific results already generated remain valid and reproducible.