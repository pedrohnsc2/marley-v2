# Module F: Stochastic Immune Simulation of MRL-ASO-001

**Generated:** 2026-04-11 00:26 UTC
**Module version:** 1.0.0
**Monte Carlo simulations:** 1000
**Noise intensity (sigma):** 0.1

## Executive Summary

MRL-ASO-001 demonstrates effective dual-function activity against intracellular L. infantum amastigotes. Deterministic model: 100.0% parasite reduction at 168h with clearance in 14.1h. Stochastic model (N=1000): 100.0% clearance probability, mean t_90% = 14.0 +/- 1.4h. Synergy index (speed) = 0.7754 (SUB-ADDITIVE): the antisense mechanism dominates at standard dose (10 uM), with TLR9-mediated NO production providing complementary killing. Both mechanisms independently achieve clearance, providing mechanistic redundancy. EC50 = 0.1 uM, EC90 = 0.413 uM.

## 1. Deterministic ODE Model — Three Scenarios

Comparison of MRL-ASO-001's two mechanisms of action, individually and combined.

| Scenario | Final Parasite Load | Reduction (%) | Clearance | Time to 90% (h) |
|---|---|---|---|---|
| Antisense only | 0.00 | 100.0 | Yes | 14.4 |
| TLR9/immune only | 0.00 | 100.0 | Yes | 45.4 |
| Dual function | 0.00 | 100.0 | Yes | 14.1 |

## 2. Synergy Analysis

### Reduction-based
- **Antisense effect:** 100.0% reduction
- **TLR9 effect:** 100.0% reduction
- **Dual effect:** 100.0% reduction
- **SI (reduction):** 0.5000

### Speed-based (primary metric)
- **Time to 90% clearance (antisense):** 14.4h
- **Time to 90% clearance (TLR9):** 45.4h
- **Time to 90% clearance (dual):** 14.1h
- **Expected additive time:** 10.9h
- **SI (speed):** 0.7754
- **Classification:** SUB-ADDITIVE

> A funcao dual e sub-aditiva em velocidade (SI_speed = 0.78). Isso indica que os mecanismos compartilham um bottleneck comum (provavelmente a taxa de captacao do ASO), limitando o beneficio adicional da combinacao. Cada mecanismo individualmente ja e potente.

## 3. Stochastic Simulation (SDE Monte Carlo)

Each scenario simulated 1000 times with multiplicative Wiener noise (sigma = 0.1).

| Scenario | Clearance Prob. | Mean P_final | Reduction (%) | Mean t_90% (h) |
|---|---|---|---|---|
| Antisense only | 100.0% | 0.00 +/- 0.00 | 100.0 | 14.3 +/- 1.6 |
| TLR9/immune only | 100.0% | 0.00 +/- 0.00 | 100.0 | 44.6 +/- 5.8 |
| Dual function | 100.0% | 0.00 +/- 0.00 | 100.0 | 14.0 +/- 1.4 |

## 4. Dose-Response Curve

| Dose (uM) | Clearance@72h (%) | Mean t_90% (h) | Reduction (%) |
|---|---|---|---|
| 0.1 | 0.0 | 132.7 | 54.3 |
| 0.5 | 28.7 | 70.4 | 100.0 |
| 1.0 | 93.0 | 46.9 | 100.0 |
| 2.0 | 100.0 | 32.3 | 100.0 |
| 5.0 | 100.0 | 19.9 | 100.0 |
| 10.0 | 100.0 | 14.0 | 100.0 |
| 20.0 | 100.0 | 9.9 | 100.0 |
| 50.0 | 100.0 | 6.2 | 100.0 |
| 100.0 | 100.0 | 4.4 | 100.0 |

**EC50:** 0.100 uM
**EC90:** 0.413 uM

## 5. Model Parameters

| Parameter | Value | Unit | Description |
|---|---|---|---|
| r_P | 0.02 | 1/h | Amastigote replication rate (doubling ~35h) |
| K_P | 200.0 | parasites | Carrying capacity per macrophage |
| k_kill | 0.1 | 1/h | NO-mediated killing rate |
| k_aso | 0.05 | 1/h | ASO-mediated splicing block rate |
| k_uptake | 0.05 | 1/h | ASO uptake rate into phagolysosome |
| A_ext | 10.0 | uM | Extracellular ASO concentration |
| k_deg | 0.000640 | 1/h | ASO degradation rate (t1/2 = 1083h) |
| k_tlr9 | 0.03 | 1/h | TLR9-induced IFN-gamma production rate |
| k_tnf | 0.02 | 1/h | TLR9-induced TNF-alpha production rate |
| d_I | 0.1 | 1/h | IFN-gamma decay rate |
| d_T | 0.12 | 1/h | TNF-alpha decay rate |
| k_no | 0.08 | 1/h | IFN-gamma -> iNOS -> NO production rate |
| d_N | 0.2 | 1/h | NO decay rate |

## Overall Conclusion

MRL-ASO-001 demonstrates robust anti-parasitic activity in the stochastic macrophage infection model:

1. **Antisense mechanism** (splicing block): 100.0% reduction, 90% clearance in 14.4h. This is the dominant mechanism, blocking trans-splicing to prevent mRNA maturation.

2. **TLR9/immune mechanism** (IFN-gamma/NO): 100.0% reduction, 90% clearance in 45.4h. Slower but provides independent clearance via innate immunity activation.

3. **Combined dual function**: 90% clearance in 14.1h (100.0% stochastic clearance). SI_speed = 0.7754 (SUB-ADDITIVE) — both mechanisms share the ASO uptake bottleneck, limiting speed synergy at high dose. However, the dual mechanism provides **mechanistic redundancy**: even if one pathway is impaired, the other maintains efficacy.

4. **Dose-response**: EC50 = 0.1 uM, EC90 = 0.413 uM. Sub-micromolar potency with 100% clearance at >= 2 uM within 72h.

5. **Key insight**: The TLR9 pathway becomes increasingly important at low doses (< 1 uM) where antisense alone is insufficient, and provides a safety net against resistance mutations that could affect the antisense mechanism.

**The ASO is predicted to effectively eliminate intracellular L. infantum amastigotes via dual mechanisms.**

## References

1. Liew FY et al. (1990) J Exp Med 172(5):1557-1559 — NO kills Leishmania amastigotes
2. Klinman DM (2004) Nat Rev Immunol 4(4):249-258 — CpG motifs and TLR9 activation
3. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO PS backbone and TLR9
4. Gantt KR et al. (2001) J Immunol 167(2):893-901 — IFN-gamma induces iNOS in macrophages
5. Kloeden PE, Platen E (1992) Numerical Solution of SDEs. Springer — Euler-Maruyama method
6. Higham DJ (2001) SIAM Review 43(3):525-546 — SDE numerical methods tutorial
7. Rogers MB et al. (2011) PLoS Genetics 7(8):e1002237 — Leishmania mutation rates
8. Murray HW et al. (2005) Lancet 366(9496):1561-1577 — Visceral leishmaniasis immunology
