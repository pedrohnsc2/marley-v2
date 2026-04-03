# Marley 🐾

<p align="center">
  <img src="docs/marley.png" alt="Marley" width="400">
</p>

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
01_fetch_genome      — Downloads all annotated protein sequences (~8,500)
        ↓
02_filter_surface     — Filters surface-exposed proteins via SignalP 6.0
        ↓
03_conservation       — Scores conservation across Brazilian strains via BLAST
        ↓
04_immunogenicity     — Predicts canine MHC binding via IEDB + loads pre-validated antigens
        ↓
05_report             — Generates ranked candidate list + Markdown report
        ↓
06_construct          — Designs multi-epitope mRNA vaccine construct
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
| 06_construct | ✅ Complete | mRNA vaccine construct designer (see below) |
| Web dashboard | ✅ MVP | Next.js + Tailwind, live data from Supabase |
| CI/CD | ✅ Complete | GitHub Actions: lint (ruff) + pytest on Python 3.11/3.12/3.13 |
| Test suite | ✅ 61 tests | Models, pipeline modules, construct designer, validated antigens |

### End-to-end validation

The full pipeline has been validated on real data:

| Stage | Input | Output | Time |
|---|---|---|---|
| Genome fetch | TriTrypDB API | 8,527 protein sequences (8.35 MB) | ~2s |
| Surface filter | 8,527 proteins | Signal peptide candidates (Sec/SPI) | ~77s/50 proteins |
| Conservation | Surface proteins | Candidates with >80% identity across strains | ~62s/protein |
| Immunogenicity | Conserved candidates | IC50 binding scores for 3 DLA alleles | ~12s/protein |
| Report | All scored candidates | Ranked list + validated antigens from literature | <1s |
| Construct design | Top candidates | mRNA vaccine sequence ready for synthesis | ~30s |

**Full analysis report:** [`docs/marley_full_report.md`](docs/marley_full_report.md) — complete results from the pipeline run on the entire *L. infantum* proteome (8,527 proteins → 15 epitopes → mRNA construct).

---

## mRNA Vaccine Construct Design (Module 06)

Module 06 takes the ranked antigen candidates and designs a complete multi-epitope mRNA vaccine construct:

### Construct architecture

```
[tPA signal peptide] → [L7/L12 adjuvant] → EAAAK → [CTL epitopes joined by AAY] → GPGPG → [HTL epitopes joined by GPGPG]
```

| Component | Description |
|-----------|-------------|
| **Signal peptide** | tPA leader (MDAMKRGLCCVLLLCGAVFVSAS) — drives secretion for MHC presentation |
| **Adjuvant** | 50S ribosomal L7/L12 — TLR4 agonist with Th1 bias (critical for *Leishmania*) |
| **CTL epitopes** | 9-mer peptides selected from IEDB MHC-I predictions (IC50 < 500 nM, 3 DLA alleles) |
| **Linkers** | AAY (proteasomal cleavage), GPGPG (prevents junctional neoepitopes), EAAAK (rigid spacer) |

### What the module produces

| Output | Description |
|--------|-------------|
| `results/construct/vaccine_construct.fasta` | Multi-epitope protein sequence |
| `results/construct/vaccine_mrna.fasta` | Full mRNA: 5'UTR + codon-optimized ORF + 3'UTR(x2) + poly(A)120 |
| `results/construct/construct_card.json` | Identity card: MW, pI, instability index, GRAVY, GC content |
| `results/construct/construct_report.md` | Design rationale with epitope table, safety assessment, references |

### Design features

- **Codon optimization** for *Canis lupus familiaris* (Kazusa codon usage table)
- **Restriction site removal** (EcoRI, BamHI, HindIII, XbaI, NheI, BsaI)
- **Homopolymer breaking** (no runs > 4 nt)
- **Physicochemical analysis** via Biopython ProtParam (MW, pI, instability, GRAVY)
- **Antigenicity check** via VaxiJen 2.0 (threshold > 0.4)
- **Allergenicity check** via AllerTOP v2.0
- **Configurable** via CLI: `--signal-peptide tPA|IgK` and `--adjuvant L7L12|RS09`

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
