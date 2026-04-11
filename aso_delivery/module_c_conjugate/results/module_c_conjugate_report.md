# Module C: Delivery Conjugate Analysis for MRL-ASO-001

**Generated:** 2026-04-11 02:20 UTC
**Module version:** 1.0.0
**ASO:** MRL-ASO-001 (25 nt, LNA-DNA-LNA gapmer, full PS)

## Executive Summary

Module C evaluated five conjugation strategies for macrophage-targeted delivery of MRL-ASO-001. **Trimannose-ASO (trivalent cluster)** is the recommended strategy with a composite score of 0.715 and expected 9.7x uptake improvement over naked ASO. The mannose receptor (MRC1/CD206) pathway delivers cargo directly to lysosomes/phagolysosomes where L. infantum amastigotes reside, eliminating the need for endosomal escape. GalNAc conjugation was confirmed UNSUITABLE: ASGPR is hepatocyte-exclusive with zero expression on macrophages.

## 1. Mannose Conjugation (MRC1/CD206 Pathway)

The macrophage mannose receptor (MRC1/CD206) is highly expressed on
macrophages and delivers cargo directly to lysosomes/phagolysosomes.

### Receptor Profile

- **Receptor:** Mannose receptor C-type 1 (MRC1 (CD206))
- **Expression on macrophages:** 8.5/10
- **Expression on hepatocytes:** 0.1/10
- **Selectivity (mac/hep):** 85.0x
- **Pathway:** Clathrin-mediated endocytosis -> early endosome -> lysosome
- **Endpoint:** lysosome/phagolysosome

### Monovalent vs Trivalent Mannose

| Property | Monovalent | Trivalent (cluster) |
|---|---|---|
| Ligand | D-mannose | alpha-1,3/alpha-1,6-trimannose |
| MW total (Da) | 8544 | 8926 |
| logP estimated | -3.472 | -3.576 |
| Kd for MRC1 (uM) | 5.0 | 0.1 |
| Uptake fold vs naked | 2.5x | 9.7x |
| Synthesis complexity | 2/5 | 4/5 |
| Cost per mg (USD) | $150 | $450 |

**Avidity improvement (tri/mono):** 50.0x in Kd, 3.9x in uptake

**Lysosomotropic advantage:**
The mannose receptor (MRC1) internalizes cargo via clathrin-mediated endocytosis, delivering directly to lysosomes. In Leishmania-infected macrophages, phagolysosomes are the compartment where amastigotes reside. This pathway places the ASO in direct proximity to the SL RNA target without requiring endosomal escape.

## 2. Cholesterol Conjugation (LDLR Pathway)

- **MW total:** 8751 Da
- **logP:** -2.899 (delta vs naked: +0.601)
- **Membrane association:** 2.00x
- **Kd for LDLR:** 2.0 uM
- **Uptake fold:** 2.5x
- **Macrophage selectivity:** poor

> **Warning:** Cholesterol conjugation primarily drives hepatic accumulation via LDLR, which is 2.8x more expressed on hepatocytes than macrophages. In vivo, cholesterol-ASO will associate with circulating LDL, directing ~60-80% of the dose to the liver. Macrophage delivery is achievable but not the primary distribution pathway.

## 3. GalNAc Assessment — NEGATIVE RESULT

**Conclusion: UNSUITABLE for macrophage delivery.**

GalNAc (N-acetylgalactosamine) is the industry gold standard for
hepatocyte-targeted ASO/siRNA delivery (Alnylam platform). However,
it is fundamentally incompatible with macrophage targeting.

### Reasons for Exclusion

1. ASGPR (asialoglycoprotein receptor) is expressed EXCLUSIVELY on hepatocytes
1. Expression on macrophages is zero — no receptor-mediated uptake possible
1. GalNAc-ASO would be directed to the liver, not to infected macrophages
1. High synthesis cost ($800/mg) with zero therapeutic benefit for this indication
1. The extraordinary affinity (Kd = 3 nM) is wasted on a non-target cell type

**Clinical context:** GalNAc is the gold standard for hepatocyte delivery (Alnylam platform: givosiran, inclisiran, lumasiran, vutrisiran). However, this platform is fundamentally incompatible with macrophage-targeted ASO delivery. The ASGPR receptor is hepatocyte-specific by design — it clears desialylated glycoproteins from circulation, a function unique to the liver.

**Note:** This negative result is documented to prevent future re-investigation. No modification of GalNAc valency, linker, or attachment site can overcome the fundamental absence of ASGPR on macrophages.

## 4. Conjugate Comparison Matrix

| Conjugate | Receptor | Expression | Kd (uM) | Uptake (fold) | Complexity | Cost/mg | Score | Suitable |
|---|---|---|---|---|---|---|---|---|
| Trimannose-ASO (trivalent cluster) | MRC1 (CD206) | 8.5/10 | 0.100 | 9.7x | 4/5 | $450 | 0.715 | Yes |
| Mannose-ASO (monovalent) | MRC1 (CD206) | 8.5/10 | 5.000 | 2.5x | 2/5 | $150 | 0.676 | Yes |
| Palmitate-ASO | SR-A (CD204) | 7.0/10 | 5.000 | 2.3x | 1/5 | $80 | 0.611 | Yes |
| Cholesterol-ASO | LDLR | 3.2/10 | 2.000 | 2.5x | 2/5 | $120 | 0.395 | **No** |
| GalNAc-ASO (trivalent) | ASGPR | 0.0/10 | 0.003 | 0.1x | 5/5 | $800 | 0.226 | **No** |
| Naked ASO (PS backbone) | Non-specific / SR-A | 7.0/10 | N/A | 1.0x | 0/5 | $50 | 0.420 | Yes |

### Ranking

1. **Trimannose-ASO (trivalent cluster)** (score: 0.715) — SUITABLE: MRC1 (CD206) is expressed on macrophages (8.5/10) with high selectivity vs hepatocytes. Strong receptor affinity (Kd = 0.100 uM). Expected uptake: 9.7x over naked ASO. Internalization via Clathrin-mediated endocytosis -> early endosome -> lysosome.
2. **Mannose-ASO (monovalent)** (score: 0.676) — SUITABLE: MRC1 (CD206) is expressed on macrophages (8.5/10) with high selectivity vs hepatocytes. Moderate receptor affinity (Kd = 5.0 uM). Expected uptake: 2.5x over naked ASO. Internalization via Clathrin-mediated endocytosis -> early endosome -> lysosome.
3. **Palmitate-ASO** (score: 0.611) — SUITABLE: MSR1 (CD204) is expressed on macrophages (7.0/10) with moderate selectivity vs hepatocytes. Moderate receptor affinity (Kd = 5.0 uM). Expected uptake: 2.3x over naked ASO. Internalization via Non-clathrin endocytosis -> endosome -> lysosome.
4. **Cholesterol-ASO** (score: 0.395) — UNSUITABLE: LDLR has poor macrophage selectivity (macrophage/hepatocyte ratio = 0.36).
5. **GalNAc-ASO (trivalent)** (score: 0.226) — UNSUITABLE: ASGPR1/ASGPR2 is not expressed on macrophages (expression = 0.0/10). This receptor is hepatocyte lysosome (liver-specific)-specific.

## 5. Optimal Conjugate Design

### Primary Recommendation

**Trimannose-ASO (trivalent cluster)**

- **Structure:** MRL-ASO-001 (25-nt LNA-DNA-LNA gapmer, full PS) conjugated at 3' end via C6 aminohexyl linker to trimannose moiety
- **Target receptor:** MRC1 (CD206)
- **Expected uptake:** 9.7x over naked ASO
- **Score:** 0.715
- **Synthesis complexity:** 4/5
- **Cost:** $450/mg

### Alternative

**Mannose-ASO (monovalent)** (score: 0.676)
- Expected uptake: 2.5x
- If Trimannose-ASO (trivalent cluster) synthesis proves too complex or costly, Mannose-ASO (monovalent) offers a simpler alternative with 2.5x uptake.

### Dual Conjugation Option

- **Primary:** Trimannose cluster (3' end)
- **Secondary:** Cy5 fluorescent tag (5' end)
- **MW total:** 9718 Da
- **Purpose:** Dual conjugation enables simultaneous macrophage targeting (trimannose at 3') and real-time tracking of ASO uptake and intracellular distribution (Cy5 at 5'). The Cy5 tag does not interfere with RNase H activity because the 5' LNA flank provides steric protection.
- **Caveat:** Dual conjugation increases MW substantially and may reduce cellular uptake. Should be used for mechanism-of-action studies, not as the therapeutic candidate.

### Final Rationale

Trimannose cluster conjugation is recommended as the primary delivery strategy for MRL-ASO-001 based on three converging advantages: (1) MRC1/CD206 is highly expressed on macrophages (8.5/10) with 85x selectivity over hepatocytes; (2) the mannose receptor internalization pathway delivers cargo directly to lysosomes/phagolysosomes, eliminating the need for endosomal escape; (3) multivalent mannose binding achieves nanomolar affinity (Kd ~100 nM) through the avidity effect, a 50-fold improvement over monovalent mannose. This strategy exploits the unique biology of Leishmania-infected macrophages: the same phagocytic receptors that the parasite exploits for entry can be co-opted for therapeutic ASO delivery.

## Overall Conclusion

MRL-ASO-001 should be conjugated with a **trimannose cluster** at the 3' end via a C6 aminohexyl linker for optimal macrophage delivery. This strategy exploits three converging biological advantages:

1. **Receptor expression:** MRC1/CD206 is highly expressed on macrophages (8.5/10) with 85x selectivity over hepatocytes, ensuring preferential macrophage uptake.

2. **Lysosomotropic delivery:** The mannose receptor pathway delivers cargo directly to the lysosome/phagolysosome compartment where L. infantum amastigotes reside. The ASO arrives at the same subcellular location as its SL RNA target.

3. **Avidity effect:** Trivalent mannose improves Kd from ~5 uM (monovalent) to ~100 nM (50-fold), achieving efficient receptor engagement at therapeutic concentrations.

The GalNAc platform (industry gold standard for hepatocyte delivery) was evaluated and definitively excluded: ASGPR is not expressed on macrophages. Cholesterol conjugation offers moderate uptake but poor macrophage selectivity due to high hepatic LDLR expression.

**Next step:** Conjugate synthesis and in vitro uptake validation in canine macrophage cell line (DH82).

## References

1. East L, Isacke CM (2002) Biochim Biophys Acta 1572(2-3):364-386 — Mannose receptor biology
2. Taylor ME et al. (1992) J Biol Chem 267(3):1719-1726 — Multivalent mannose binding
3. Martinez-Pomares L (2012) J Leukoc Biol 92(6):1177-1186 — MRC1 on macrophages
4. Wolfrum C et al. (2007) Nat Biotechnol 25(10):1149-1157 — Cholesterol-siRNA delivery
5. Nair JK et al. (2014) JACS 136(49):16958-16961 — GalNAc-siRNA platform
6. Spiess M (1990) Biochemistry 29(43):10009-10018 — ASGPR liver specificity
7. Toulme JJ et al. (1994) Biochimie 76(3-4):153-154 — Palmitate-ASO macrophages
8. Kawakami S et al. (2008) J Control Release 131(3):153-158 — Mannosylated liposomes Leishmania
9. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review
10. Goldstein JL, Brown MS (2009) Arterioscler Thromb Vasc Biol 29(4):431-438 — LDL receptor
11. Peiser L et al. (2002) Curr Opin Immunol 14(1):123-128 — Scavenger receptor macrophages
12. Irache JM et al. (2008) J Control Release 128(1):15-25 — Mannosylated nanoparticles
