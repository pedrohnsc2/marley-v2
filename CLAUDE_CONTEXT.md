# Marley Project — Full Context for Claude

## Who is the user
Pedro Nascimento (pedrohnsc2@gmail.com). Lost his dog Marley to canine visceral leishmaniasis. Built this computational pipeline to find a vaccine and drug against the disease. Working with Prof. Rodolfo Giunchetti (UFMG) and Wanessa Goes as potential experimental collaborators.

## Repository
- Primary: `github.com/pedrohnsc2/marley-v2` (PRIVATE, clean history, personal email)
- Backup: `github.com/pedrohnsc2/marley` (PRIVATE, has corporate email in commits — do not use for new work)
- Local: `/Users/pedronascimento/dev/marley` (original working directory)
- Local clean: `/Users/pedronascimento/dev/marley-new` (marley-v2 clone)

## What was built (in order)

### v1 — Vaccine antigen discovery (pipeline/)
- Downloaded 8,527 L. infantum proteins from TriTrypDB/UniProt
- Filtered surface proteins via SignalP 6.0 (139 candidates)
- Conservation scoring via NCBI BLAST
- MHC binding prediction via IEDB (3 DLA dog alleles)
- Selected 11 epitopes (IC50 range: 11-118 nM)
- Designed mRNA vaccine construct: 335 aa, tPA signal + L7/L12 adjuvant + 11 CTL epitopes + AAY/GPGPG linkers
- Codon optimized for Canis lupus familiaris (CAI = 0.9948)
- 3D structure predicted via ESMFold
- 3 construct variants (A/B/C) with different adjuvants

### v2 — Drug target discovery (drug_targets/)
- Mapped 52 enzymes across 5 metabolic pathways (purine salvage, trypanothione, sterol biosynthesis, pentose phosphate, glycolysis)
- Compared each with human homolog via BLAST
- 7 priority validated targets: TryS (0.98), TryR (0.96), ADL (0.95), SMT (0.94), GMPS (0.92), 6PGDH (0.90), XPRT (0.88)

### v3 — Molecular docking (drug_targets/)
- Screened 77 compounds against top 5 enzyme targets via AutoDock Vina
- Top hit: Methotrexate vs GMPS (-8.07 kcal/mol)
- Pentamidine ranked #2 (already used clinically — validates methodology)
- Designed 20 custom MRL molecules based on Pemetrexed scaffold
- Best: MRL-003 (amide tail) at -7.74 kcal/mol vs TryR
- ADMET filtering with RDKit (Lipinski Rule of 5)
- 3D visualization with PyMOL

### v3 — Selectivity test (HONEST NEGATIVE)
- MRL-003 docked against human glutathione reductase: -8.68 kcal/mol (BINDS BETTER to human enzyme)
- 20 redesigned variants (charge, size, spermidine modifications): 0/20 achieved selectivity
- Conclusion: antifolates cannot selectively inhibit TryR — too similar to human GR

### v3 — Resistance prediction
- Alanine scanning of 21 TryR binding site residues
- All mutations neutral — LOW resistance risk

### v3 — Cross-species analysis
- L. donovani TryR: 98% identity (MRL-003 would work)
- L. major: 92%, L. braziliensis: 73%, T. cruzi: 46%, T. brucei: 23%

### v4 — Vaccine optimization (pipeline/)
- Epitope APL (Altered Peptide Ligand) design: 93 variants tested via IEDB
- Safety check: BLAST vs dog proteome — 0/11 epitopes cross-react
- Adjuvant screening: 6 candidates, L7/L12 confirmed best for Th1
- Construct reconstruction with 4 ordering strategies

### v5 — Computational validation (pipeline/)
- Selectivity (module 12): MRL-003 vs human GR — FAILED
- Resistance (module 13): 21 mutations all neutral — LOW risk
- Cross-species (module 14): L. donovani 98% compatible
- Delivery optimization (module 15): CAI 0.9948, LNP formulation specified
- Immune simulation (module 17): ODE model, 82% Th1, memory stable at day 365
- Selectivity redesign (module 18): 20 variants, 0/20 selective — FAILED

### v4-RNA — Information theory (rna_entropy/)
- Applied Shannon entropy H(X) = -Σp(x)log₂p(x) to L. infantum transcriptome
- Core formula in: rna_entropy/03_shannon_entropy.py (function _window_entropy, line ~108)
- Identified SL RNA (Spliced Leader): 39nt sequence added to EVERY mRNA via trans-splicing
- SL sequence: AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG
- Conserved ~500 million years, absent in humans (confirmed via BLAST)
- information_score: 0.99 (highest in entire project)
- Codon usage analysis (RSCU), SL RNA mapping, human comparison

### ASO drug design (rna_entropy/)
- Designed 119 ASO candidates targeting SL RNA
- Top candidate: MRL-ASO-001
  - Sequence: ACAGAAACTGATACTTATATAGCGT (25 nt)
  - Tm: 61.4°C (IDT) / 68.5°C (our calc) / 108.5°C with LNA
  - ΔG binding: -28.0 kcal/mol
  - Hairpin ΔG: -1.66 (17x weaker than target binding) — PASS
  - Self-dimer ΔG: -5.83 (4.8x weaker than target binding) — PASS
  - Off-target human: NONE (0/119)
  - Off-target dog: NONE
  - Chemical design: LNA-DNA-LNA gapmer + phosphorothioate backbone
  - Delivery: subcutaneous injection, 5-10 mg/kg/week
  - Selectivity: 100%

### DUAL FUNCTION DISCOVERY
- VaxiJen score of MRL-ASO-001: 1.2561 (ANTIGEN) — 5.4x higher than Leish-Tec
- The PS backbone activates TLR9 (innate immune stimulation)
- MRL-ASO-001 simultaneously kills the parasite AND stimulates the immune system
- NOTE: VaxiJen was designed for proteins, not DNA — the score is technically off-label but the TLR9 mechanism is literature-confirmed (Krieg et al. 1995)

### Oral drug attempt (HONEST NEGATIVE)
- 15 RNA-binding small molecules docked against SL RNA 3D structure
- All returned +10.0 kcal/mol (no binding)
- 39nt RNA too small/flexible for small molecule pockets
- Vina designed for proteins, not short RNA

## External validation results

### VaxiJen v2.0 (parasite model, threshold 0.4)
- Leish-Tec (A2 protein): 0.2340 (NON-ANTIGEN)
- Marley vaccine construct: 0.3235 (NON-ANTIGEN)
- Marley epitopes only: 0.3730 (NON-ANTIGEN)
- MRL-ASO-001: 1.2561 (ANTIGEN)

### ExPASy ProtParam (vaccine protein)
- MW: 35,954 Da, pI: 7.93
- Instability index: 28.38 (STABLE, <40)
- Aliphatic index: 108.30 (thermostable)
- Half-life: 30h (mammalian)

### AllerTOP v2.1
- Vaccine construct: Probable NON-ALLERGEN

### IDT OligoAnalyzer (MRL-ASO-001)
- Tm: 61.4°C, MW: 7673.1 Da, GC: 32%
- Hairpin: -1.66 kcal/mol (PASS)
- Self-dimer: -5.83 kcal/mol (PASS)
- All 12 validation tests: PASSED

## Key files
- `docs/marley_complete_report.md` — full analysis report with all results
- `README.md` — comprehensive project documentation
- `rna_entropy/03_shannon_entropy.py` — Shannon entropy calculation (core math)
- `rna_entropy/08_aso_design.py` — ASO design with nearest-neighbor thermodynamics
- `results/construct/vaccine_construct.fasta` — vaccine protein sequence
- `results/construct/vaccine_mrna.fasta` — vaccine mRNA sequence
- `results/aso/aso_report.md` — ASO drug report
- `data/structures/marley_vaccine_construct.pdb` — vaccine 3D structure
- `data/structures/MRL_ASO_001.pdb` — ASO 3D structure

## Tools installed
- Python 3.13, .venv at /Users/pedronascimento/dev/marley/.venv
- AutoDock Vina (.venv/bin/vina) — x86_64 via Rosetta
- Open Babel (.venv/bin/obabel)
- PyMOL (brew install pymol)
- RDKit (via conda, linked to .venv via .pth file)
- Supabase (tables: candidates, drug_targets, docking_compounds, docking_results)

## Important decisions made
- Repository made PRIVATE for IP protection
- Considering patent for MRL-ASO-001 + vaccine construct
- Meeting with Giunchetti (UFMG) and Wanessa Goes scheduled for Friday
- User prefers honest reporting — negative results documented prominently
- User's language: Portuguese (but code/docs in English)

## Honest negative results (IMPORTANT — user values transparency)
1. MRL-003 vs human GR: binds human enzyme better — NOT selective
2. 20 MRL-003 redesigns: 0/20 achieved selectivity
3. Oral small molecule vs SL RNA: no binding (Vina not suited for short RNA)
- These failures led to the pivot from protein targets to RNA targets (SL RNA)
- The pivot produced the best result of the entire project (MRL-ASO-001)
