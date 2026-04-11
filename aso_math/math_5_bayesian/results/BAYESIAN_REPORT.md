# Bayesian Optimization of ASO Design Parameters

**Generated:** 2026-04-11 01:31 UTC
**Runtime:** 1.77 seconds
**Total evaluations:** 251
**Protocol:** 50 initial + 200 BO iterations

---

## Executive Summary

**Bayesian optimization found a potentially better variant.** The best design (2-14-14, composite improvement: 0.073641) outperforms MRL-ASO-001 in the composite objective. However, this result requires experimental validation — the composite score is a simplified mathematical model.

MRL-ASO-001 ranks **#197** out of 251 designs evaluated (percentile: 21.6%).

## 1. Design Space

| Parameter | Range | Description |
|-----------|-------|-------------|
| ASO length | 20-30 nt | Total oligonucleotide length |
| Gap size | 10-20 nt | Central DNA gap for RNase H |
| Start position | 0-14 | Position on SL RNA (for 25-mer) |
| LNA asymmetry | [0, 1] | Distribution of LNA flanks |

**Constraint:** LNA modifications only in flanking regions (gapmer design).
**Minimum LNA flank:** 2 nt per side.

## 2. Multi-Objective Function

| Objective | Weight | Goal | Source |
|-----------|--------|------|--------|
| dG binding | 35% | Minimize (stronger binding) | NN model |
| Hairpin dG | 15% | Maximize (less self-structure) | NN model |
| Self-dimer dG | 15% | Maximize (less dimerization) | NN model |
| Tm penalty | 20% | Zero in 55-75 C range | NN model |
| Off-target | 15% | Minimize (fewer off-targets) | GC/length proxy |

## 3. Gaussian Process Surrogate Model

- **Kernel:** RBF (Squared Exponential)
- **Final length scale:** 0.3000
- **Final signal variance:** 2.0000
- **Noise variance:** 1.00e-04
- **Log marginal likelihood:** 35.75
- **Hyperparameter optimization:** Grid search on log marginal likelihood
- **Refit frequency:** Every 25 iterations

## 4. Top 5 Variants Found

| Rank | Sequence | Arch | Start | Composite | dG (kcal/mol) | Tm (C) | Hairpin dG | Self-dimer dG | Is MRL? |
|------|----------|------|-------|-----------|---------------|--------|------------|--------------|---------|
| 1 | `ACAGAAACTGATACTTATATAGCGTTAGTT` | 2-14-14 | 0 | -0.186691 | -33.34 | 71.3 | 0.00 | -2.79 | no |
| 2 | `ACAGAAACTGATACTTATATAGCGTTAGTT` | 8-20-2 | 0 | -0.186691 | -33.34 | 71.3 | 0.00 | -2.79 | no |
| 3 | `ACAGAAACTGATACTTATATAGCGTTAGTT` | 10-18-2 | 0 | -0.186691 | -33.34 | 71.3 | 0.00 | -2.79 | no |
| 4 | `ACAGAAACTGATACTTATATAGCGTTAGTT` | 6-13-11 | 0 | -0.186691 | -33.34 | 71.3 | 0.00 | -2.79 | no |
| 5 | `ACAGAAACTGATACTTATATAGCGTTAGTT` | 10-16-4 | 0 | -0.186691 | -33.34 | 71.3 | 0.00 | -2.79 | no |

## 5. Detailed Comparison: Best Found vs MRL-ASO-001

| Metric | Best Found | MRL-ASO-001 | Difference |
|--------|-----------|-------------|------------|
| Sequence | `ACAGAAACTGATACTTATATAGCGTTAGTT` | `ACAGAAACTGATACTTATATAGCGT` | — |
| Architecture | 2-14-14 | 4-17-4 | — |
| Composite score | -0.186691 | -0.113050 | 0.073641 |
| dG binding | -33.34 | -27.97 | 5.37 kcal/mol |
| Hairpin dG | 0.00 | 0.00 | 0.00 kcal/mol |
| Self-dimer dG | -2.79 | -2.56 | 0.23 kcal/mol |
| Tm | 71.3 C | 68.5 C | -2.9 C |

## 6. Pareto Optimality Assessment

- **Pareto-optimal designs (5D):** 144 out of 251
- **MRL-ASO-001 is Pareto-optimal:** NO
- **Distance to Pareto front:** 0.224400 (normalized objective space)

### Pareto Front: dG Binding vs Hairpin Penalty (2D projection)

| Arch | f1 (binding) | f2 (hairpin) | Composite | Is MRL? |
|------|-------------|-------------|-----------|---------|
| 6-13-11 | -0.8335 | -0.0000 | -0.186691 | no |
| 10-16-4 | -0.8335 | -0.0000 | -0.186691 | no |
| 10-14-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 3-10-17 | -0.8335 | -0.0000 | -0.177841 | no |
| 6-16-8 | -0.8335 | -0.0000 | -0.186691 | no |
| 12-12-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 9-18-3 | -0.8335 | -0.0000 | -0.186691 | no |
| 8-20-2 | -0.8335 | -0.0000 | -0.186691 | no |
| 13-12-5 | -0.8335 | -0.0000 | -0.177841 | no |
| 15-11-4 | -0.8335 | -0.0000 | -0.186691 | no |
| 7-16-7 | -0.8335 | -0.0000 | -0.177841 | no |
| 6-19-5 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-20-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 6-18-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 2-11-17 | -0.8335 | -0.0000 | -0.177841 | no |
| 9-16-5 | -0.8335 | -0.0000 | -0.186691 | no |
| 13-12-5 | -0.8335 | -0.0000 | -0.186691 | no |
| 12-10-8 | -0.8335 | -0.0000 | -0.186691 | no |
| 6-19-5 | -0.8335 | -0.0000 | -0.177841 | no |
| 3-17-10 | -0.8335 | -0.0000 | -0.186691 | no |
| 6-17-7 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-18-8 | -0.8335 | -0.0000 | -0.177841 | no |
| 16-10-4 | -0.8335 | -0.0000 | -0.186691 | no |
| 14-11-5 | -0.8335 | -0.0000 | -0.186691 | no |
| 5-19-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 2-19-9 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-20-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-18-8 | -0.8335 | -0.0000 | -0.186691 | no |
| 9-10-11 | -0.8335 | -0.0000 | -0.177841 | no |
| 6-13-11 | -0.8335 | -0.0000 | -0.177841 | no |
| 2-18-10 | -0.8335 | -0.0000 | -0.186691 | no |
| 9-12-9 | -0.8335 | -0.0000 | -0.177841 | no |
| 2-19-9 | -0.8335 | -0.0000 | -0.177841 | no |
| 10-18-2 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-20-6 | -0.8335 | -0.0000 | -0.177841 | no |
| 8-20-2 | -0.8335 | -0.0000 | -0.186691 | no |
| 15-12-3 | -0.8335 | -0.0000 | -0.186691 | no |
| 17-10-3 | -0.8335 | -0.0000 | -0.177841 | no |
| 12-12-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 3-19-8 | -0.8335 | -0.0000 | -0.186691 | no |
| 2-11-17 | -0.8335 | -0.0000 | -0.186691 | no |
| 14-10-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 5-15-10 | -0.8335 | -0.0000 | -0.186691 | no |
| 12-16-2 | -0.8335 | -0.0000 | -0.186691 | no |
| 4-20-6 | -0.8335 | -0.0000 | -0.186691 | no |
| 9-19-2 | -0.8335 | -0.0000 | -0.186691 | no |
| 3-13-14 | -0.8335 | -0.0000 | -0.186691 | no |
| 5-14-11 | -0.8335 | -0.0000 | -0.186691 | no |
| 2-15-13 | -0.8335 | -0.0000 | -0.186691 | no |
| 8-20-2 | -0.8335 | -0.0000 | -0.177841 | no |
| 2-14-14 | -0.8335 | -0.0000 | -0.186691 | no |

## 7. Convergence Diagnostics

- **Initial best composite:** -0.175195
- **Final best composite:** -0.186691
- **Total improvement:** 0.011496
- **Converged at iteration:** ~9

## 8. Limitations & Caveats

1. **Simplified model:** The composite objective uses nearest-neighbor thermodynamics and heuristic off-target scoring. Real ASO efficacy depends on cellular uptake, RNase H recruitment, metabolic stability, and other factors not modeled here.
2. **LNA thermodynamics not modeled:** LNA modifications increase Tm by ~3-5 C per modification, but the NN model uses DNA/DNA parameters. This means absolute Tm values are underestimates for LNA gapmers.
3. **4D parameter space:** The search explores a continuous relaxation of the discrete design space. Some promising designs in intermediate parameter regions may have been missed.
4. **No experimental validation:** All results are computational predictions. Any promising variant must be synthesized and tested in vitro before drawing conclusions about superiority.

## 9. Methods

**Surrogate model:** Gaussian Process with RBF kernel, implemented from scratch.
**Acquisition function:** Expected Improvement (EI) with xi=0.01.
**Hyperparameter optimization:** Grid search maximizing log marginal likelihood.
**Seed designs:** 50 random + MRL-ASO-001 reference.
**Optimization iterations:** 200.
**EI maximization:** Random sampling with 5000 candidates per iteration.
**Pareto analysis:** Exhaustive pairwise dominance check (5 objectives).

**References:**
- Rasmussen & Williams (2006) *Gaussian Processes for Machine Learning*
- Jones, Schonlau & Welch (1998) *Efficient Global Optimization*
- SantaLucia (1998) *PNAS* 95(4):1460-1465
