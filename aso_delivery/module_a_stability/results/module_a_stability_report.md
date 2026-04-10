# Module A: Phagolysosomal Stability Analysis of MRL-ASO-001

**Generated:** 2026-04-10 23:50 UTC
**Module version:** 1.0.0

## Executive Summary

MRL-ASO-001 PASSES all phagolysosomal stability tests. At pH 4.5: dG = -23.97 kcal/mol (threshold: -15.0 kcal/mol), fraction bound = 1.0000, gapmer half-life = 1083.04 hours, LNA geometry maintained at 99.0% C3'-endo. The LNA-DNA-LNA gapmer with full PS backbone is well-suited for the hostile phagolysosomal environment of L. infantum amastigotes.

## 1. Thermodynamic Stability Across pH Gradient

The ASO traverses a pH gradient from blood (7.4) to phagolysosome (4.5).
Protonation of cytosine N3 (pKa 4.2) and adenine N1 (pKa 3.5)
destabilizes Watson-Crick base pairs at low pH.

| Compartment | pH | dG (kcal/mol) | Correction | Tm (C) | Fraction Bound | Functional |
|---|---|---|---|---|---|---|
| blood_plasma | 7.4 | -27.99 | +0.0069 | 108.5 | 1.0000 | Yes |
| early_endosome | 6.5 | -27.95 | +0.0549 | 108.3 | 1.0000 | Yes |
| late_endosome | 5.0 | -26.45 | +1.5543 | 103.8 | 1.0000 | Yes |
| phagolysosome | 4.5 | -23.97 | +4.0345 | 96.4 | 1.0000 | Yes |

**Verdict:** MRL-ASO-001 **PASSES** the pH stability test.
dG at pH 4.5 remains well below the functional threshold of -15.0 kcal/mol.

## 2. Nuclease Resistance at pH 4.5

First-order degradation kinetics: [ASO](t) = [ASO]_0 * exp(-k*t)

| Backbone | k (1/h) | Half-life (h) | f(24h) | f(48h) | Therapeutic Window |
|---|---|---|---|---|---|
| PO (phosphodiester) | 0.500000 | 1.39 | 0.000006 | 0.000000 | Not met |
| PS (phosphorothioate) | 0.001000 | 693.15 | 0.976286 | 0.953134 | Met |
| LNA-DNA-LNA gapmer (PS backbone) | 0.000640 | 1083.04 | 0.984757 | 0.969747 | Met |

**PS/PO resistance ratio:** 500.0x
**Gapmer/PO resistance ratio:** 781.2x

**Gapmer architecture:** 5-15-5 (LNA-DNA-LNA)
**Effective k:** 0.00064000/h
**Estimated half-life:** 1083.04 hours

**Verdict:** MRL-ASO-001 **PASSES** the nuclease resistance test.

## 3. LNA Conformational Stability

LNA nucleotides maintain C3'-endo sugar pucker via covalent methylene bridge,
independent of pH. This is critical for high-affinity hybridization.

| pH | Compartment | LNA C3'-endo | DNA C3'-endo | Geometry | Penalty (kcal/mol) |
|---|---|---|---|---|---|
| 7.4 | blood_plasma | 99.0% | 36.0% | Maintained | 0.0000 |
| 6.5 | early_endosome | 99.0% | 33.3% | Maintained | 0.2025 |
| 5.0 | late_endosome | 99.0% | 28.8% | Maintained | 0.5400 |
| 4.5 | phagolysosome | 99.0% | 27.3% | Maintained | 0.6525 |

**LNA Tm boost:** 4.71 C per LNA residue (10 residues total)
**Verdict:** MRL-ASO-001 **PASSES** the LNA conformational stability test.

## 4. Comparison with FDA-Approved ASOs

| Property | Nusinersen | Mipomersen | Inotersen | MRL-ASO-001 |
|---|---|---|---|---|
| Length (nt) | 18 | 20 | 20 | 25 |
| Backbone | full PS | full PS | full PS | full PS |
| Gapmer design | N/A (uniform 2'-MOE, splice-switching) | 5-10-5 | 5-10-5 | 5-15-5 |
| Route | intrathecal | subcutaneous | subcutaneous | to be determined (parenteral expected) |
| Target pH | 7.3 | 4.5 | 4.5 | 4.5 |
| Half-life (h) | 3240.0 | 720.0 | 624.0 | 1083.04 |
| Mechanism | splice-switching (steric block) | RNase H1-mediated mRNA degradation | RNase H1-mediated mRNA degradation | RNase H1-mediated SL RNA degradation |

**Key insight:** MRL-ASO-001 has a unique advantage over subcutaneous ASOs like
mipomersen and inotersen: the target RNA (SL RNA) and the ASO colocalize in the
same acidic compartment (phagolysosome), eliminating the need for endosomal escape.

## Overall Conclusion

MRL-ASO-001 demonstrates robust stability across all four analyses:

1. **pH stability:** dG decreases from -28.0 kcal/mol (pH 7.4) to -23.97 kcal/mol (pH 4.5), a loss of only 4.03 kcal/mol. The duplex remains far below the functional threshold (-15.0 kcal/mol).

2. **Nuclease resistance:** The LNA-DNA-LNA gapmer with PS backbone has an estimated half-life of 1083.04 hours at pH 4.5, far exceeding the 24-hour therapeutic window requirement.

3. **LNA conformation:** The methylene bridge maintains C3'-endo sugar pucker at 99.0% even at pH 4.5, ensuring optimal hybridization geometry.

4. **FDA comparison:** MRL-ASO-001's chemistry (LNA gapmer, full PS) matches or exceeds the stability profiles of clinically approved ASOs like mipomersen that operate in similar acidic compartments.

**The ASO is predicted to remain stable and functional in the phagolysosomal environment of canine macrophages.**

## References

1. SantaLucia J Jr (1998) PNAS 95(4):1460-1465 — Nearest-neighbor parameters
2. Siegfried NA et al. (2010) Biochemistry 49(15):3225-3236 — Base protonation at low pH
3. Eckstein F (2014) Nucleic Acids Res 42(6):3777-3788 — Phosphorothioate stability
4. Vester B, Wengel J (2004) Biochemistry 43(42):13233-13241 — LNA conformational locking
5. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review
6. Kurreck J (2003) Eur J Biochem 270(8):1628-1644 — ASO design principles
7. Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — ASO pharmacokinetics
8. McTigue PM et al. (2004) Biochemistry 43(18):5388-5405 — LNA Tm enhancement
