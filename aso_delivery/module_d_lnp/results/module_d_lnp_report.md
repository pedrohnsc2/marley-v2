# Module D: LNP Formulation Analysis for MRL-ASO-001

**Generated:** 2026-04-11 00:23 UTC
**Module version:** 1.0.0

## Executive Summary

LNP encapsulation is RECOMMENDED for MRL-ASO-001 delivery to macrophages. Key metrics: encapsulation efficiency = 99% (N/P 8), phagolysosomal release = 100%, particle diameter = 87 nm (phagocytable), mannose-targeted uptake = 3.5x vs untargeted. The pH-responsive release of DLin-MC3-DMA (pKa 6.44) delivers the ASO cargo directly to the phagolysosome where L. infantum amastigotes and their SL RNA reside. Lyophilization enables deployment in endemic regions of Brazil without cold chain infrastructure.

## 1. LNP Composition

Standard 4-component formulation based on Onpattro (patisiran).

| Component | Name | Mol% | MW (g/mol) | Role |
|---|---|---|---|---|
| Ionizable Lipid | DLin-MC3-DMA | 50.0% | 642.1 | pH-responsive encapsulation and release |
| Helper Lipid | DSPC | 10.0% | 790.1 | Bilayer structural integrity |
| Sterol | Cholesterol | 38.5% | 386.6 | Membrane stability and fusogenicity |
| Peg Lipid | DMG-PEG2000 | 1.5% | 2693.0 | Stealth layer (prevents opsonization and RES clearance) |

**Weighted average MW:** 589.32 g/mol
**Charge at pH 7.4:** 0.049407 (near-neutral, stealth)
**Charge at pH 4.5:** 0.494324 (cationic, triggers release)
**pKa (DLin-MC3-DMA):** 6.44

## 2. N/P Ratio Optimization

ASO: 25 nt, 24 PS charges (negative).
Optimal N/P range for PS-ASOs: 4-8 (vs 6-12 for mRNA).

| N/P | EE (%) | Diameter (nm) | Zeta (mV) | Phagocytable | Optimal |
|---|---|---|---|---|---|
| 1 | 42.3 | 67.5 | -23.3 | Yes | No |
| 2 | 66.7 | 70.0 | -17.6 | Yes | No |
| 3 | 80.8 | 72.5 | -10.0 | Yes | No |
| 4 | 88.9 | 75.0 | -2.4 | Yes | **Yes** |
| 5 | 93.6 | 77.5 | 3.3 | Yes | **Yes** |
| 6 | 96.3 | 80.0 | 6.7 | Yes | **Yes** |
| 7 | 97.9 | 82.5 | 8.4 | Yes | **Yes** |
| 8 | 98.8 | 85.0 | 9.3 | Yes | **Yes** |
| 9 | 99.0 | 87.5 | 9.7 | Yes | No |
| 10 | 99.0 | 90.0 | 9.8 | Yes | No |
| 11 | 99.0 | 92.5 | 9.9 | Yes | No |
| 12 | 99.0 | 95.0 | 10.0 | Yes | No |
| 13 | 99.0 | 97.5 | 10.0 | Yes | No |
| 14 | 99.0 | 100.0 | 10.0 | Yes | No |
| 15 | 99.0 | 102.5 | 10.0 | Yes | No |
| 16 | 99.0 | 105.0 | 10.0 | Yes | No |
| 17 | 99.0 | 107.5 | 10.0 | Yes | No |
| 18 | 99.0 | 110.0 | 10.0 | Yes | No |
| 19 | 99.0 | 112.5 | 10.0 | Yes | No |
| 20 | 99.0 | 115.0 | 10.0 | Yes | No |

**Selected N/P:** 8
**Encapsulation efficiency:** 98.8%
**Particle diameter:** 85.0 nm
**Zeta potential:** 9.3 mV

## 3. pH-Responsive Release Kinetics

Model: f_released(pH) = 1 / (1 + exp(4 * (pH - 6.44)))

| Compartment | pH | MC3 Protonated (%) | ASO Released (%) | Status | Transit (min) |
|---|---|---|---|---|---|
| blood_plasma | 7.4 | 9.9 | 2.1 | minimal | 0 |
| early_endosome | 6.5 | 46.6 | 44.0 | partial | 5 |
| late_endosome | 5.0 | 96.5 | 99.7 | complete | 15 |
| phagolysosome | 4.5 | 98.9 | 100.0 | complete | 30 |

**Release at target (phagolysosome):** 100.0%

**Key insight:** LNP pH-responsive release is ADVANTAGEOUS for macrophage delivery: at pH 4.5 (phagolysosome), 100.0% of the ASO cargo is released. Since the target RNA (SL RNA) is located in the SAME compartment where L. infantum amastigotes reside, endosomal escape is NOT required — the LNP releases MRL-ASO-001 directly where it needs to act.

## 4. Macrophage-Targeted LNP (Mannose-PEG)

**Target receptor:** CD206 (mannose receptor / MRC1)

L. infantum-infected macrophages upregulate CD206 expression. Mannose-decorated LNPs are recognized by CD206 and internalized via receptor-mediated endocytosis, directing the ASO payload to the phagolysosome where the parasite resides.

| Mannose-PEG (%) | Uptake Fold | Diameter (nm) | Phagocytable |
|---|---|---|---|
| 0 | 2.0x | 85.0 | Yes |
| 25 | 2.8x | 86.0 | Yes |
| 50 | 3.5x | 87.0 | Yes |
| 75 | 4.2x | 88.0 | Yes |
| 100 | 5.0x | 89.0 | Yes |

**Recommended:** 50% mannose-PEG substitution
**Expected uptake increase:** 3.5x

## 5. Storage Stability

Arrhenius model: k(T) = k_ref * exp(Ea/R * (1/T_ref - 1/T))

| Temp (C) | k_agg (1/day) | Half-life (days) | f(30d) | f(90d) | f(365d) |
|---|---|---|---|---|---|
| -20 | 0.000203 | 3414.9 | 99.4% | 98.2% | 92.9% |
| +4 | 0.002396 | 289.2 | 93.1% | 80.6% | 41.7% |
| +25 | 0.015000 | 46.2 | 63.8% | 25.9% | 0.4% |

**25 C shelf life adequate (>80% at 90 days):** No

### Lyophilization for Tropical Deployment

- **Feasible:** Yes
- **Cryoprotectants:** sucrose (10% w/v), trehalose (10% w/v)
- **Reconstitution:** 5 min in sterile water, vortex 30 sec
- **Size increase post-reconstitution:** <20%
- **EE post-reconstitution:** >80%
- **Brazil advantage:** Eliminates cold chain dependency for deployment in Northeast Brazil (Bahia, Ceara, Maranhao) where visceral leishmaniasis is endemic and cold chain infrastructure is limited.

## 6. LNP vs Naked ASO Comparison

| Property | Naked ASO | LNP-ASO |
|---|---|---|
| Cellular uptake (relative) | 1.0x | 3.5x |
| Encapsulation efficiency | N/A | 99% |
| Target delivery | passive (gymnosis or electroporation) | mannose receptor-mediated endocytosis |
| Nuclease protection | PS backbone only | LNP encapsulation + PS backbone |
| Cost per dose (USD) | $5.00 | $15.00 |
| Storage | Lyophilized powder, room temperature | 4 C (liquid) or lyophilized |
| Manufacturing complexity | 2/5 | 4/5 |

**Effective dose reduction:** 3.5x
**Cost ratio (LNP/naked):** 3x
**Cost-effectiveness ratio:** 1.15

**Recommendation:** LNP formulation is RECOMMENDED for MRL-ASO-001 because: (1) 3.5x higher macrophage uptake via mannose-CD206, (2) 100% release at phagolysosomal pH 4.5, (3) effective dose reduction of 3.5x justifies the 3x cost increase. Cost-effectiveness ratio: 1.15.

## Overall Conclusion

LNP encapsulation provides significant advantages for MRL-ASO-001 delivery:

1. **Encapsulation:** 99% efficiency at N/P 8, protecting the ASO during circulation.

2. **pH-responsive release:** 100% of cargo released at phagolysosomal pH 4.5. Unlike mRNA/siRNA therapeutics that require endosomal escape, MRL-ASO-001 benefits from phagolysosomal release because the target (SL RNA) is in the SAME compartment.

3. **Macrophage targeting:** Mannose-PEG-DSPE decoration provides 3.5x higher uptake by infected macrophages via CD206 receptor-mediated endocytosis.

4. **Size:** 87 nm diameter is within the optimal range for phagocytic uptake (50-200 nm).

5. **Storage:** Lyophilization with sucrose/trehalose enables room-temperature storage, critical for deployment in endemic areas of Northeast Brazil.

6. **Cost-effectiveness:** Despite 3x higher cost per dose, the 3.5x uptake improvement and 99% encapsulation yield a cost-effectiveness ratio of 1.15, justifying the LNP formulation.

**The LNP-formulated MRL-ASO-001 is predicted to be a viable delivery strategy for targeting L. infantum in canine macrophages.**

## References

1. Cullis PR, Hope MJ (2017) Mol Ther 25(7):1467-1475 — LNP design principles
2. Jayaraman M et al. (2012) Angew Chem 51(34):8529-8533 — DLin-MC3-DMA optimization
3. Semple SC et al. (2010) Nat Biotechnol 28(2):172-176 — LNP encapsulation
4. Patel S et al. (2020) J Control Release 327:146-160 — Macrophage targeting with mannose
5. Kulkarni JA et al. (2018) Nanoscale 10(10):4567-4573 — pH-responsive release mechanism
6. Sahay G et al. (2013) Nat Biotechnol 31(7):653-658 — Endosomal escape and trafficking
7. Schoenmaker L et al. (2021) Int J Pharm 601:120586 — LNP storage stability
8. Ball RL et al. (2017) Drug Deliv Transl Res 7(1):89-103 — LNP lyophilization
9. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review
