# MARLEY PROJECT — COMPLETE CONTEXT FOR CLAUDE

## Who

Pedro Nascimento, software developer and independent bioinformatics researcher in Minas Gerais, Brazil. Lost dog Marley to canine visceral leishmaniasis. Built this computational drug discovery platform in his memory. Applying for PPGIT master's at UFMG (enrollment April 13, 2026).

**Collaborators:**
- Prof. Dr. Rodolfo Cordeiro Giunchetti (ICB/UFMG) — experimental collaborator, potential advisor
- Dra. Wanessa M. Goes (CTVacinas/UFMG) — bioinformatics co-collaborator
- Dra. Silvane Maria Fonseca Murta (Fiocruz Rene Rachou) — experimental validation partner

## Repository

- GitHub: `github.com/pedrohnsc2/marley` (PRIVATE)
- Local: `/Users/pedronascimento/dev/marley`
- Python 3.13, venv at `./venv`
- PyTorch 2.11 (MPS/Apple Silicon), ESM-2, Transformers, Anthropic SDK

## Project Scale

- **84,000+ lines of Python**
- **455 automated tests** (all passing)
- **80+ JSON result files** (4.2 MB total)
- **2 validation certificates** (math 52/60 + delivery 60/60)
- **3 parallel tracks** (ASO therapy + vaccine platforms + AI/ML)
- **11 functional AI modules**

---

## DISCOVERY PIPELINE (v1-v4)

### v1 — Vaccine antigen discovery (`pipeline/`, 14.5k lines)
- 8,527 L. infantum proteins from TriTrypDB/UniProt
- SignalP 6.0 → 139 surface proteins
- BLAST conservation + IEDB MHC binding (3 DLA alleles)
- **11 CTL epitopes selected** (IC50: 11-118 nM)
- mRNA construct: 335 aa, tPA signal + L7/L12 adjuvant + AAY/GPGPG linkers
- CAI = 0.9948 for Canis lupus familiaris
- 3D structure predicted via ESMFold

### v2 — Drug targets (`drug_targets/`, 6.1k lines)
- 52 enzymes in 5 metabolic pathways vs human homologs via BLAST
- 7 validated targets: TryS (0.98), TryR (0.96), ADL (0.95), SMT (0.94), GMPS (0.92), 6PGDH (0.90), XPRT (0.88)

### v3 — Molecular docking
- 77 compounds vs 5 targets via AutoDock Vina
- Pentamidine #2 → methodology validated
- 20 custom MRL molecules (Pemetrexed scaffold)
- **HONEST NEGATIVE: MRL-003 binds human GR better (-8.68 vs -7.74 kcal/mol). 0/20 redesigns selective. Pivoted to RNA.**

### v4-RNA — Information theory (`rna_entropy/`, 8.3k lines)
- Shannon entropy H(X) on full L. infantum transcriptome
- **SL RNA identified** (information_score = 0.99)
- 39nt, trans-spliced onto every mRNA, conserved ~500M years, absent in mammals
- 119 ASO candidates → **MRL-ASO-001 selected**
- Zero off-targets in human or dog (BLAST confirmed)
- Dual function: antisense (RNase H) + TLR9 immunostimulation

---

## TRACK 1 — MRL-ASO-001 (Therapeutic)

### MRL-ASO-001 Profile
- Type: 25-nt LNA-DNA-LNA gapmer, phosphorothioate backbone
- Target: L. infantum Spliced Leader RNA (positions 5-30)
- Sequence: ACAGAAACTGATACTTATATAGCGT (CONFIDENTIAL — never expose externally)
- dG binding: -27.97 kcal/mol
- Tm: 68.48 C
- Off-targets: 0 (human + canine)
- Dual mechanism: antisense + TLR9

### Mathematical Validation (`aso_math/`, 16.7k lines) — 52/60 VALIDATED
Certificate: `aso_math/reports/results/math_certificate_v2.json`

1. Thermodynamic Optimality: 6/10 — dG -27.97, rank 45/76 mutants, optimal at 30nt but 25nt is delivery trade-off
2. Information Geometry (Fisher): 10/10 — 14.5x host distance ratio, alienness 0.1607
3. Topological Stability (TDA): 10/10 — 5 persistent H1 features, topology gain 5.0
4. Target Irreplaceability: 10/10 — SL RNA is #1/17 nodes, 62.6% lambda2 drop on removal
5. Design Optimization (Bayesian): 6/10 — NOT Pareto-optimal (rank 197/251), justified by length trade-off
6. Resistance Barrier (Markov): 10/10 — 0/75 escape mutations, 285 years worst-case

### Delivery Feasibility (`aso_delivery/`, 11.9k lines) — 60/60 VALIDATED
Certificate: `aso_delivery/reports/results/delivery_certificate.json`

A. Phagolysosomal Stability: 10/10 — 100% bound pH 4.5, half-life 1,083h, LNA C3'-endo 99%
B. Membrane Permeation: 10/10 — 200,246 nM intracellular (2002x threshold), macrophage advantage 18.7x
C. Conjugate Delivery: 10/10 — Trimannose-ASO via MRC1, uptake 9.7x
D. LNP Formulation: 10/10 — 98.8% encapsulation, 87nm, 100% release pH 4.5
E. ADMET Profile: 10/10 — TI 8.0, bioavailability 87%, t1/2 21 days
F. Immune SDE: 10/10 — 100% clearance (N=1000), EC50 0.1 uM, clearance 14h

### Integration Report
`results/integration/aso_integration_report.json` + `ASO_INTEGRATION_REPORT.md`
Cross-track synthesis of math + delivery + therapy comparison with 7 documented limitations.

---

## TRACK 2 — Vaccine Platforms (`vaccine_platforms/`, 9.5k lines)

3 platforms from the same 11 immutable epitopes (`shared/epitopes.py`):

| Platform | Cost/dose (industrial) | vs Leish-Tec | Time to market |
|----------|----------------------|-------------|----------------|
| A — mRNA-LNP | $5.75 | -85.8% | 36-60 months |
| B — E. coli (pET-28a) | $2.28 | -91.6% | 24-36 months |
| C — L. tarentolae (BSL-1) | $0.00008 | -100% | 31-48 months |

**Critical finding:** Platform C SL RNA is identical to L. infantum in ASO binding region (dG = 0). Combined ASO + live vaccine requires temporal separation >2 weeks.

Market: R$880.75M/year (54.2M dogs in Brazil, 35.2M at risk, 14.1M doses/year)

Comparison report: `vaccine_platforms/reports/results/comparison_results.json`

---

## TRACK 3 — AI/ML (`marley_ai/`, 8.7k lines)

11 functional modules. Config: frozen dataclasses. Results: envelope JSON. Device: MPS (Apple M3 Pro, 36GB).

| Module | Status | Key Result |
|--------|--------|-----------|
| 01_rag | Working | 288 PubMed papers, TF-IDF + Claude Haiku for synthesis |
| 02_leish_kg | Working | 43 nodes, 58 edges (NetworkX) |
| 03_leish_esm | Working | ESM-2 650M, 25 sequences, 1280-dim embeddings, epitopes cluster together |
| 04_rna_fm | Working | Custom RNA transformer, Nussinov structure, 13/13 pairs disrupted by ASO, 94.9% conservation |
| 05_rosettafold | Working | 3D duplex, dG -27.38 (2.1% from experimental), RNase H accessible, 3 PDB files |
| 06_evodiff | Working | Discrete diffusion, 32 ASO + 47 epitope variants, best dG -29.25 |
| 07_contrastive | Working | InfoNCE, 100% allele accuracy, Spearman -0.19 vs IC50 |
| 08_rl_ppo | Working | REINFORCE, +28.5% reward, best dG -36.83 (but GC 60% = off-target risk) |
| 09_sae | Working | 10 interpretable features, hydrophobicity dominates (7/10), top selectivity 46x |
| 10_digital_twin | Working | Coupled PK+ODE+SDE, 28-day simulation of infected dog |
| 11_scientist | Working | 5 agents, 5 hypotheses, 5 experiments, Claude synthesis, consensus 0.60 |

### AI Scientist Top Hypotheses
1. (85%) MRL-ASO-001 inhibits trans-splicing via RNase H
2. (75%) RL variants stronger binding but off-target risk at GC > 55%
3. (70%) N-terminal hydrophobicity predicts antigenicity in Leishmania
4. (65%) Optimal ASO design: GC 40-50%, dG -30 to -33
5. (60%) LINF_230010300 (hypothetical protein) is virulence factor

### AI Scientist Top Experiments
1. RT-qPCR trans-splicing in promastigotes (wet lab, priority 1)
2. Evaluate RL variants with RNA-FM + RoseTTAFold (computational)
3. Generate GC 0.40-0.50 ASO variants with EvoDiff (computational)
4. ELISPOT canino with top 5 epitopes (wet lab)
5. Annotate LINF_230010300 via structural homology (computational)

---

## TESTS

455 tests, all passing in <15s:
- `tests/aso_delivery/` — 202 tests (modules A-F: stability, membrane, conjugate, LNP, ADMET, SDE)
- `tests/vaccine_platforms/` — 99 tests (shared epitopes, platforms A/B/C)
- `tests/aso_math/` — 66 tests (integration, modules, organisms, target_config, thermo)
- `tests/` — 48 tests (pipeline, core, construct, scoring, structure)

Run: `./venv/bin/python -m pytest tests/ -v`

---

## KEY FILES

| File | Description |
|------|------------|
| `aso_math/config.py` | SSOT: ASO_SEQUENCE, SL_SEQUENCE, NN_PARAMS |
| `aso_math/thermo.py` | compute_dg(), compute_tm(), gc_content() |
| `aso_math/reports/results/math_certificate_v2.json` | Math validation 52/60 |
| `aso_delivery/reports/results/delivery_certificate.json` | Delivery validation 60/60 |
| `aso_delivery/run_all.py` | Orchestrator for 6 delivery modules |
| `vaccine_platforms/shared/epitopes.py` | 11 immutable epitopes (SSOT) |
| `vaccine_platforms/reports/comparison_matrix.py` | 3-platform comparison |
| `marley_ai/config.py` | AIModuleConfig, detect_device(), MODULES registry |
| `marley_ai/registry.py` | @register decorator, AIModule Protocol |
| `marley_ai/11_scientist/discovery_loop.py` | Discovery loop + Claude synthesis |
| `results/integration/aso_integration_report.json` | Cross-track synthesis |
| `results/construct/construct_card.json` | Vaccine construct identity card |
| `docs/MEDIUM_ARTICLE.md` | Medium article draft |
| `README.md` | Complete project documentation |

---

## ENVIRONMENT

- `.env` contains: SUPABASE_URL, SUPABASE_KEY, ANTHROPIC_API_KEY (all gitignored)
- Claude API: connected (Haiku for RAG/Scientist, ~$0.003/call)
- Hardware: Apple M3 Pro, 36GB RAM, MPS available
- PyTorch 2.11, transformers 5.5.3, fair-esm 2.0.0, anthropic 0.94.0

---

## RULES

1. **11 epitopes are IMMUTABLE** (`vaccine_platforms/shared/epitopes.py`)
2. **verify_sl_rna_orthogonality()** mandatory before Platform C work
3. **NEVER expose externally:** MRL-ASO-001 sequence or vaccine_mrna.fasta
4. **Document negative results** — they're data, not failures
5. **Comments in Portuguese**, code/variables in English
6. **All results in JSON** in `results/` for reproducibility
7. **rna_entropy/03_shannon_entropy.py** and **08_aso_design.py** — DO NOT MODIFY
8. Module 11 (Multi-Agent) only after all other modules working

---

## HONEST NEGATIVE RESULTS

1. MRL-003 binds human GR better than parasite TryR (-8.68 vs -7.74). 0/20 redesigns selective.
2. 15 oral compounds vs SL RNA 3D: none bound. ASO is only viable approach.
3. MRL-ASO-001 NOT Pareto-optimal (rank 197/251). 25nt is deliberate trade-off.
4. Dual mechanism synergy is sub-additive (0.78), not synergistic as hypothesized.
5. Platform C SL RNA identical to L. infantum in binding region. Can't combine ASO + live vaccine.
6. VaxiJen 2.0 classified Leish-Tec as NON-ANTIGEN (0.2340) despite 96% clinical efficacy. Tool limitation.
7. RL top variants all have GC 60% — off-target risk flagged by AI Scientist.

---

## WHAT'S NEXT

All results are computational. Zero wet-lab data. Priority experimental validation:
1. Synthesize MRL-ASO-001 (~$500)
2. RT-qPCR in promastigotes (validates mechanism)
3. ITC/SPR binding (confirms dG -27.97)
4. DH82 macrophage assay (validates intracellular killing)
5. Platform B preclinical in BALB/c
6. Off-target transcriptomics (RNA-seq)
7. Canine PK study (validates 87% bioavailability)
