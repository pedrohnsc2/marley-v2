# Marley

> In memory of Marley, lost to visceral leishmaniasis.
> This pipeline is dedicated to every dog that didn't make it.

An open-source bioinformatics pipeline for **reverse vaccinology** targeting canine visceral leishmaniasis (*Leishmania infantum*).

---

## What it does

Marley identifies mRNA vaccine antigen candidates through a five-stage computational pipeline. Starting from the complete *L. infantum* proteome, it progressively filters, scores, and ranks proteins to surface the most promising targets for a canine mRNA vaccine -- no wet-lab work required until the final candidates are selected.

The goal is to lower the barrier to vaccine research for a disease that kills thousands of dogs every year and remains endemic in Latin America, the Mediterranean, and parts of Asia.

---

## Pipeline stages

| Stage | Name | Description |
|-------|------|-------------|
| 1 | **Genome Fetch** | Downloads the *L. infantum* proteome from [TriTrypDB](https://tritrypdb.org). Caches the FASTA locally so subsequent runs skip the download. |
| 2 | **Surface Filtering** | Identifies surface-exposed proteins using SignalP 6.0 signal-peptide prediction. Only proteins likely to be accessible to the host immune system pass this gate. |
| 3 | **Conservation Analysis** | Evaluates protein conservation across *Leishmania* strains via NCBI BLAST. Highly conserved proteins make better vaccine targets because they are less likely to mutate and escape immune recognition. |
| 4 | **Immunogenicity Scoring** | Predicts MHC-I and MHC-II binding affinities through the IEDB Analysis Resource API. Proteins that generate strong T-cell responses score higher. |
| 5 | **Report Generation** | Combines all scores into a weighted final ranking, persists results to Supabase, and generates a human-readable report of the top antigen candidates. |

---

## Stack

| Component | Version / Notes |
|-----------|-----------------|
| Python | 3.11+ |
| Biopython | Sequence parsing and BLAST interaction |
| requests / httpx | HTTP clients for TriTrypDB, IEDB, and SignalP APIs |
| pandas | Tabular data manipulation and report generation |
| tqdm | Progress bars for long-running stages |
| python-dotenv | Environment variable management |
| supabase-py | Supabase client for persisting candidates |
| pytest | Test framework |

---

## Setup

```bash
# 1. Clone the repository
git clone https://github.com/your-username/marley.git
cd marley

# 2. Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate   # Linux / macOS
# .venv\Scripts\activate    # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure environment variables
cp .env.example .env
# Edit .env and fill in your Supabase URL and anon key:
#   SUPABASE_URL=https://your-project.supabase.co
#   SUPABASE_KEY=your-anon-key
```

### Supabase setup

1. Create a free project at [supabase.com](https://supabase.com).
2. Open the SQL Editor and run the schema in [`docs/supabase_schema.sql`](docs/supabase_schema.sql).
3. Copy the project URL and anon key into your `.env` file.

---

## How to run

```bash
# Run the full pipeline (all 5 stages)
python run_pipeline.py

# Skip the genome download (use cached proteome)
python run_pipeline.py --skip-fetch

# Dry run -- execute all stages but do not write to the database
python run_pipeline.py --dry-run

# Combine flags
python run_pipeline.py --skip-fetch --dry-run
```

### Output

- A ranked table of antigen candidates printed to stdout.
- Results persisted to the `candidates` table in Supabase (unless `--dry-run` is set).
- A CSV report written to `output/candidates.csv`.

---

## Database schema (Supabase)

The pipeline stores its results in a single `candidates` table. Full DDL is in [`docs/supabase_schema.sql`](docs/supabase_schema.sql).

```sql
create table candidates (
  id                    uuid default gen_random_uuid() primary key,
  gene_id               text unique not null,
  gene_name             text,
  sequence              text,
  has_signal_peptide    boolean default false,
  conservation_score    float,
  immunogenicity_score  float,
  final_score           float,
  filters_passed        text[],
  status                text default 'pending',
  created_at            timestamp default now(),
  updated_at            timestamp default now()
);
```

Key indexes:

- `idx_candidates_gene_id` -- fast lookup by gene identifier.
- `idx_candidates_final_score` -- efficient ranking queries (descending).

---

## Contributing

Contributions are welcome. To get started:

1. Fork the repository and create a feature branch from `main`.
2. Write tests for any new functionality (`pytest`).
3. Make sure all tests pass before opening a pull request.
4. Follow the existing code style (the project uses standard Python conventions).
5. Open a pull request with a clear description of what you changed and why.

If you find a bug or have a feature request, please open an issue first so we can discuss the approach.

---

## Scientific context

### Reverse vaccinology

Traditional vaccine development starts in the lab -- growing pathogens, extracting proteins, and testing them one by one. **Reverse vaccinology** flips this process: it starts with the pathogen's genome and uses computational tools to predict which proteins are most likely to provoke a protective immune response. Only the best candidates are then synthesized and tested experimentally.

This approach was pioneered by Rino Rappuoli for *Neisseria meningitidis* serogroup B and has since been applied to dozens of pathogens.

### Visceral leishmaniasis

Visceral leishmaniasis (VL), also known as kala-azar, is a vector-borne disease caused by *Leishmania* parasites transmitted through sandfly bites. In dogs, the disease is often fatal and serves as the primary urban reservoir for human infection. There is currently no widely available, highly effective vaccine for canine VL.

### Why mRNA vaccines for dogs

mRNA vaccines offer several advantages for veterinary use:

- **Speed**: once antigens are identified, mRNA constructs can be synthesized in weeks.
- **Safety**: no live pathogen is involved; the mRNA degrades naturally.
- **Flexibility**: multi-antigen constructs can target several proteins simultaneously.
- **Scalability**: manufacturing is simpler than traditional protein-based vaccines.

By combining reverse vaccinology with mRNA technology, Marley aims to accelerate the development of an effective canine VL vaccine.

---

## License

MIT
