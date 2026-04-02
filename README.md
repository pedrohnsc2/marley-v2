# Marley 🐾

> *In memory of Marley, lost to canine visceral leishmaniasis.*
> *This pipeline is dedicated to every dog that didn't make it.*

[![Status](https://img.shields.io/badge/status-in%20development-yellow)]()
[![License](https://img.shields.io/badge/license-MIT-blue)]()
[![Python](https://img.shields.io/badge/python-3.11+-green)]()

An open-source bioinformatics pipeline for reverse vaccinology targeting canine visceral leishmaniasis (*Leishmania infantum*). Marley automates the computational identification and ranking of antigen candidates for mRNA vaccine development.

---

## What it does

Marley connects to public genomic databases and runs a multi-stage filtering pipeline to identify which proteins from *L. infantum* are the strongest candidates for a canine mRNA vaccine — combining computational prediction with a curated list of experimentally validated antigens from Brazilian researchers.

---

## Pipeline stages

```
TriTrypDB (L. infantum genome)
        ↓
01_fetch_genome      — Downloads all annotated protein sequences (~8,000)
        ↓
02_filter_surface     — Filters surface-exposed proteins via SignalP 6.0
        ↓
03_conservation       — Scores conservation across Brazilian strains via BLAST
        ↓
04_immunogenicity     — Predicts canine MHC binding via IEDB + loads pre-validated antigens
        ↓
05_report             — Generates ranked candidate list + Markdown report
```

---

## Experimentally validated antigens

The following antigens are pre-loaded with priority status, sourced from published research by Brazilian groups (UFMG, UFOP, Fiocruz/MG):

| Antigen | Source | Evidence | Score |
|---|---|---|---|
| LiHyp1 | Giunchetti/UFMG | Murine validation, Th1 response via IFN-γ, immunoproteomics | 0.95 |
| A2 | UFMG / Leish-Tec | Only MAPA-approved vaccine in Brazil, 96.41% efficacy | 0.92 |
| LBSap_antigens | Reis/UFOP + Giunchetti/UFMG | Technology transferred to Ouro Fino Saúde Animal | 0.90 |
| Lutzomyia_longipalpis_proteins | Giunchetti/UFMG | Patent052, transmission-blocking | 0.88 |
| KMP-11 | Literature | High conservation across strains, documented immunogenicity | 0.85 |
| LiESP_Q | Literature | High diagnostic specificity, immunoprotective potential | 0.83 |
| LACK | Literature | T-cell activation in murine models | 0.82 |
| HSP70_HSP83 | Literature | Overexpressed under stress, conserved, strong cellular response | 0.80 |

---

## Stack

| Layer | Technology |
|---|---|
| Core language | Python 3.11+ |
| Bioinformatics | Biopython |
| Data | pandas |
| Database | Supabase |
| Web dashboard | Next.js + TypeScript |
| CI/CD | GitHub Actions |
| Signal peptide | SignalP 6.0 via BioLib SDK |
| External APIs | TriTrypDB, NCBI BLAST, IEDB |

---

## Setup

```bash
git clone https://github.com/pedrohnsc2/marley
cd marley
pip install -r requirements.txt
cp .env.example .env
# Fill in SUPABASE_URL and SUPABASE_KEY
```

Run the Supabase schema:
```sql
create table candidates (
  id uuid default gen_random_uuid() primary key,
  gene_id text unique not null,
  gene_name text,
  sequence text,
  has_signal_peptide boolean default false,
  conservation_score float,
  immunogenicity_score float,
  final_score float,
  filters_passed text[],
  status text default 'pending',
  priority boolean default false,
  source text,
  evidence text,
  created_at timestamp default now(),
  updated_at timestamp default now()
);

create index idx_candidates_final_score on candidates(final_score desc);
create index idx_candidates_priority on candidates(priority desc);
```

---

## How to run

```bash
# Full pipeline (interactive, confirms each stage)
python run_pipeline.py

# Full pipeline (automated, saves progress between stages)
python run_full_pipeline.py

# Skip genome download if already fetched
python run_pipeline.py --skip-fetch

# Dry run (no external API calls)
python run_pipeline.py --dry-run

# Run tests
python -m pytest tests/ -v
```

### Web dashboard

```bash
cd web
npm install
cp .env.local.example .env.local
# Fill in Supabase credentials
npm run dev
# Open http://localhost:3000
```

---

## Module status

| Module | Status | Notes |
|---|---|---|
| 01_fetch_genome | ✅ Complete | 8,527 proteins downloaded from TriTrypDB |
| 02_filter_surface | ✅ Complete | SignalP 6.0 via BioLib SDK, tested end-to-end |
| 03_conservation | ✅ Complete | NCBI BLAST against L. donovani, L. major, L. braziliensis |
| 04_immunogenicity | ✅ Complete | IEDB MHC-I with 3 canine DLA alleles (netmhcpan_ba) |
| 05_report | ✅ Complete | Ranked CSV + Markdown report with validated antigens |
| Web dashboard | ✅ MVP | Next.js + Tailwind, live data from Supabase |
| CI/CD | ✅ Complete | GitHub Actions: lint (ruff) + pytest on Python 3.11/3.12/3.13 |
| Test suite | ✅ 40 tests | Models, pipeline modules, validated antigens, logger |

### End-to-end validation

The full pipeline has been validated on real data:

| Stage | Input | Output | Time |
|---|---|---|---|
| Genome fetch | TriTrypDB API | 8,527 protein sequences (8.35 MB) | ~2s |
| Surface filter | 8,527 proteins | Signal peptide candidates (Sec/SPI) | ~77s/50 proteins |
| Conservation | Surface proteins | Candidates with >80% identity across strains | ~62s/protein |
| Immunogenicity | Conserved candidates | IC50 binding scores for 3 DLA alleles | ~12s/protein |
| Report | All scored candidates | Ranked list + validated antigens from literature | <1s |

---

## Scientific context

Canine visceral leishmaniasis affects millions of dogs in Brazil, particularly in Minas Gerais. The only approved vaccine (Leish-Tec, developed at UFMG) was suspended in 2023 due to quality deviations in the A2 protein antigen. No replacement is currently available.

mRNA vaccines represent a promising path forward — the same platform that enabled rapid COVID-19 vaccine development could be applied to *Leishmania* using well-characterized antigens and lipid nanoparticle delivery. The main challenge remains inducing cellular (Th1) rather than humoral (Th2) immunity.

Marley contributes to this effort by automating the computational layer of antigen discovery, making the pipeline reproducible and open to collaboration with wet-lab researchers.

---

## Contributing

Marley is open to collaboration from both developers and researchers.

- **Developers:** check open issues, improve pipeline modules, build the web dashboard
- **Researchers:** validate computational outputs, suggest antigens, share experimental data

---

## License

MIT — free to use, modify and distribute with attribution.
