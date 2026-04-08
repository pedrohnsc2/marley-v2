# Module 06: mRNA Vaccine Construct Designer — Implementation Plan

**Status:** Approved for implementation
**Date:** 2026-04-02

---

## Overview

After modules 01-05 produce a ranked list of antigen candidates, module 06 translates those results into a tangible mRNA vaccine construct design — ready for synthesis and experimental validation.

---

## Architecture

### New files

| File | Purpose |
|------|---------|
| `pipeline/06_construct.py` | Main module (~500-600 lines) |
| `core/codon_tables.py` | Canis lupus familiaris codon usage data (Kazusa, taxid 9615) |
| `core/models.py` (extend) | Add `Epitope` and `VaccineConstruct` dataclasses |
| `core/db.py` (extend) | Add `upsert_construct` / `get_construct` helpers |
| `docs/supabase_schema_v2.sql` | New tables: `vaccine_constructs`, `construct_epitopes` |
| `tests/test_construct.py` | Unit tests |

### No new dependencies

Biopython (`Bio.SeqUtils.ProtParam`) covers physicochemical calculations. `requests` covers web API calls. `beautifulsoup4` for VaxiJen/AllerTOP HTML parsing (add to requirements.txt).

---

## Pipeline Flow

```
Top candidates (from module 05)
        ↓
1. Select best epitopes
   ├── CTL: 9-mer from IEDB MHC-I (already have)
   └── HTL: 15-mer from IEDB MHC-II (new API call)
        ↓
2. Assemble multi-epitope protein
   Signal peptide (tPA) + Adjuvant (L7/L12) + EAAAK
   + [CTL1-AAY-CTL2-AAY-...-CTLn]
   + GPGPG
   + [HTL1-GPGPG-HTL2-GPGPG-...-HTLn]
        ↓
3. Codon optimization (Canis lupus familiaris)
   ├── Top codon per amino acid
   ├── GC content 40-60%
   ├── Remove homopolymers > 4
   └── Remove restriction sites
        ↓
4. mRNA assembly
   5'UTR + ATG + optimized ORF + STOP + 3'UTR + 3'UTR + poly(A)120
        ↓
5. Safety checks
   ├── VaxiJen antigenicity (> 0.4 = desirable)
   ├── AllerTOP allergenicity (non-allergen = desirable)
   └── ProtParam: MW, pI, instability, GRAVY
        ↓
6. Output
   ├── vaccine_construct.fasta (protein)
   ├── vaccine_mrna.fasta (mRNA)
   ├── construct_card.json (identity card)
   └── construct_report.md (design rationale)
```

---

## Task Breakdown

### Phase 1 — Foundation (can run in parallel)

| Task | Complexity | Agent |
|------|-----------|-------|
| 1. Add `Epitope` + `VaccineConstruct` dataclasses to `core/models.py` | S | Implementer |
| 2. Add DB helpers for construct tables in `core/db.py` | S | Implementer |
| 3. Create `core/codon_tables.py` with canine codon frequencies | S | Implementer |
| 4. Create Supabase schema for new tables | S | Implementer |

### Phase 2 — Core logic (parallel after Phase 1)

| Task | Complexity | Agent |
|------|-----------|-------|
| 5. Epitope selection algorithm (greedy, non-overlapping, diverse) | M | Implementer |
| 6. MHC-II prediction via IEDB API (HTL epitopes) | M | Implementer |
| 7. Construct assembly (linkers, signal peptide, adjuvant) | S | Implementer |
| 8. Codon optimization + reverse translation | M | Implementer |
| 9. mRNA assembly (UTRs, poly-A) | S | Implementer |

### Phase 3 — Integration (sequential after Phase 2)

| Task | Complexity | Agent |
|------|-----------|-------|
| 10. Physicochemical analysis (ProtParam) | S | Implementer |
| 11. VaxiJen + AllerTOP integration (web scraping with fallback) | M | Implementer |
| 12. Output file generation (FASTA, JSON, Markdown) | M | Implementer |
| 13. Wire into `run_full_pipeline.py` | S | Implementer |

### Phase 4 — Testing + polish

| Task | Complexity | Agent |
|------|-----------|-------|
| 14. Write comprehensive tests | M | Tester |
| 15. End-to-end validation | M | QA |

---

## Key Design Decisions

### Construct structure

```
[tPA signal peptide] - [L7/L12 adjuvant] - EAAAK - [CTL epitopes joined by AAY] - GPGPG - [HTL epitopes joined by GPGPG]
```

| Component | Sequence | Purpose |
|-----------|----------|---------|
| tPA leader | MDAMKRGLCCVLLLCGAVFVSAS | Secretion signal |
| L7/L12 adjuvant | MAKLSTDELLD...TVTVK (120 aa) | TLR4 agonist, Th1 bias |
| EAAAK | EAAAK | Rigid linker, adjuvant → epitopes |
| AAY | AAY | CTL linker, proteasomal cleavage |
| GPGPG | GPGPG | HTL linker, prevents junctional neoepitopes |

### Epitope selection criteria

- CTL (9-mer): IC50 < 500 nM from MHC-I predictions (3 DLA alleles)
- HTL (15-mer): IC50 < 500 nM from MHC-II predictions (DLA class II or HLA-DRB1 proxy)
- Non-overlapping: skip peptides with >= 5 shared residues from same protein
- Diverse: max 3 epitopes per source protein
- Target: 12-15 CTL + 8-10 HTL epitopes

### mRNA structure

```
5'cap (m7G, annotated) + 5'UTR (alpha-globin + Kozak) + ORF + TGA + 3'UTR (beta-globin x2) + poly(A)120
```

---

## Supabase Schema Changes

```sql
CREATE TABLE vaccine_constructs (
    construct_id        TEXT PRIMARY KEY,
    protein_sequence    TEXT NOT NULL,
    mrna_sequence       TEXT NOT NULL,
    signal_peptide      TEXT NOT NULL,
    adjuvant_name       TEXT NOT NULL,
    epitope_count       INTEGER NOT NULL,
    gc_content          REAL,
    molecular_weight    REAL,
    isoelectric_point   REAL,
    instability_index   REAL,
    gravy               REAL,
    vaxijen_score       REAL,
    allergenicity       TEXT,
    created_at          TIMESTAMPTZ DEFAULT now()
);

CREATE TABLE construct_epitopes (
    id                  SERIAL PRIMARY KEY,
    construct_id        TEXT REFERENCES vaccine_constructs(construct_id),
    peptide             TEXT NOT NULL,
    epitope_type        TEXT NOT NULL,
    source_gene_id      TEXT,
    allele              TEXT NOT NULL,
    ic50                REAL NOT NULL,
    position_in_construct INTEGER
);
```

---

## Risk Register

| Risk | Severity | Mitigation |
|------|----------|------------|
| VaxiJen/AllerTOP have no REST API (web-form only) | High | Web scraping with graceful fallback; non-blocking |
| IEDB MHC-II may not support canine DLA class II | High | Fall back to HLA-DRB1*01:01 as cross-species proxy |
| Validated antigens have empty sequences | High | Fetch from UniProt (A2=Q25327, KMP-11=P15711) or skip |
| Junctional neoepitopes from linker boundaries | Medium | Post-assembly scan of junction 9-mers against IEDB |
| Codon optimization GC drift | Low | Post-processing pass with synonymous codon swaps |

---

## MVP vs Nice-to-haves

### MVP (ship first)
- CTL epitope selection from MHC-I data (already have)
- Multi-epitope construct with AAY linkers, tPA signal, L7/L12 adjuvant
- Codon optimization (top codon per AA + restriction site removal)
- Full mRNA sequence
- ProtParam identity card
- VaxiJen antigenicity (with fallback)
- Design report in Markdown

### V2 (later)
- HTL epitopes via MHC-II predictions
- AllerTOP allergenicity
- Multiple adjuvant variants for comparison
- CAI (Codon Adaptation Index) scoring
- RNAfold mRNA structure prediction
- 3D structure via AlphaFold/ESMFold
- B-cell linear epitope prediction (BepiPred)

---

## Open Questions

1. **B-cell epitopes**: Leave BCL block empty for MVP, add BepiPred later?
2. **Validated antigen sequences**: Fetch from UniProt or rely only on pipeline candidates?
3. **Multiple construct variants**: Generate one default or compare tPA vs IgK / L7L12 vs RS09?
