# Module B: Membrane Permeability and Macrophage Uptake of MRL-ASO-001

**Generated:** 2026-04-11 02:20 UTC
**Module version:** 1.0.0

## Executive Summary

MRL-ASO-001 PASSES the membrane permeability assessment. While passive diffusion is thermodynamically impossible (membrane barrier >> 20 kT), macrophage uptake via receptor-mediated endocytosis (scavenger receptors) achieves an estimated intracellular concentration of 200246 nM in the phagolysosome at 1 uM extracellular dose. This exceeds the effective threshold of 100 nM. Macrophages provide a 19x uptake advantage over non-phagocytic cells due to high SR-A expression and constitutive pinocytosis. Critically, no endosomal escape is needed because the target (SL RNA) resides in the same phagolysosomal compartment where the ASO accumulates after endocytic uptake.

## 1. Physicochemical Properties

ASOs are large (>8 kDa), highly charged (charge ~-24) macromolecules.
They violate all Lipinski criteria for passive membrane permeation.

| Property | PO (unmodified) | PS (phosphorothioate) | MRL-ASO-001 (LNA gapmer PS) |
|---|---|---|---|
| MW (Da) | 8634.0 | 8634.0 | 8774.0 |
| logP | -24.00 | -18.00 | -17.80 |
| logD (pH 7.4) | -36.00 | -30.00 | -29.80 |
| logD (pH 4.5) | -35.99 | -29.99 | -29.79 |
| Net charge (pH 7.4) | -24.0 | -24.0 | -24.0 |
| Net charge (pH 4.5) | -24.0 | -24.0 | -24.0 |
| Hydrodynamic radius (nm) | 2.72 | 2.72 | 2.72 |
| Lipinski violations | 4 | 4 | 4 |

**Conclusion:** MRL-ASO-001 has MW 8774 Da, net charge -24 at pH 7.4, and logD -29.8. This places it firmly outside the drug-likeness space for passive membrane permeation (4/4 Lipinski violations). Active uptake mechanisms are required.

## 2. Membrane Interaction Energetics

Three energy terms govern ASO-membrane interaction:
- **Electrostatic repulsion:** negative ASO vs negative membrane surface
- **Hydrophobic insertion:** PS backbone has partial hydrophobic character
- **Born solvation:** cost of transferring charges from water to lipid

| Energy term | PO | PS | MRL-ASO-001 |
|---|---|---|---|
| Electrostatic (kcal/mol) | 0.52 | 0.52 | 0.52 |
| Hydrophobic (kcal/mol) | 7.20 | -2.40 | -2.40 |
| Born solvation (kcal/mol) | 171.01 | 171.01 | 171.01 |
| Total barrier (kcal/mol) | 178.73 | 169.13 | 169.13 |
| Barrier (kT) | 290.00 | 274.40 | 274.40 |

**Partition coefficient (log K):** PO = -125.9, PS = -119.2, MRL-ASO-001 = -119.2

**Conclusion:** The total membrane barrier for MRL-ASO-001 is 169.1 kcal/mol (274 kT). Passive diffusion is infeasible (threshold: 20 kT). The Born solvation term dominates, reflecting the enormous cost of transferring 24 negative charges from water (dielectric 74) to lipid (dielectric 2). Active endocytosis is the ONLY viable uptake mechanism.

## 3. Endocytosis Pathways in Macrophages

Macrophages are professional phagocytes with 5-10x higher endocytic activity
than most cell types. Three uptake pathways are modeled:

| Pathway | Rate (1/h) | Efficiency (24h) | [Intracellular] (nM) | Dominant? |
|---|---|---|---|---|
| Fluid-phase pinocytosis | 0.005000 | 11.3% | 1615.4 | No |
| Receptor-mediated endocytosis (scavenger receptors) | 0.030000 | 51.3% | 197015.2 | Yes |
| Gymnosis (free uptake) | 0.005000 | 11.3% | 1615.4 | No |

**Total efficiency (24h):** 61.7%
**Total intracellular concentration:** 200246.1 nM
**Effective threshold:** 100 nM
**Verdict:** MRL-ASO-001 **PASSES** the intracellular concentration test.
**Macrophage advantage:** 18.7x vs non-phagocytic cells

## 4. Comparison with Clinically Approved ASOs (Delivery)

| Property | Nusinersen | Mipomersen | Inotersen | MRL-ASO-001 |
|---|---|---|---|---|
| Route | intrathecal | subcutaneous | subcutaneous | to be determined (parenteral expected) |
| Target cell | motor neurons (CNS) | hepatocytes | hepatocytes | macrophages (L. infantum-infected) |
| Uptake mechanism | Gymnosis + receptor-mediated in CNS cells | ASGPR-mediated endocytosis (liver-specific) + scavenger r... | Similar to mipomersen: ASGPR + scavenger receptors. Enhan... | Receptor-mediated endocytosis via SR-A (dominant). Gymnos... |
| Endosomal escape needed | Yes | Yes | Yes | **No** |
| Formulation | Preservative-free solution in artificial CSF | Aqueous solution for SC injection | Aqueous solution for SC injection | TBD (naked ASO may be sufficient for macrophage delivery) |

### Key Advantages of MRL-ASO-001 for Macrophage Delivery

- PS backbone provides natural tropism for macrophages via SR-A
- No endosomal escape needed (target is in phagolysosome)
- Macrophages are professional phagocytes (5-10x higher uptake)
- Infected macrophages may have enhanced endocytic activity
- LNA flanks provide nuclease resistance in acidic phagolysosome

## Overall Conclusion

MRL-ASO-001 demonstrates a favorable delivery profile for macrophage-targeted therapy against L. infantum:

1. **Physicochemical profile:** MW 8774 Da, charge -24 at pH 7.4 confirms that passive membrane permeation is impossible, as expected for ASOs.

2. **Membrane barrier:** The total energy barrier exceeds 20 kT, dominated by Born solvation energy. This is consistent with the known requirement for active uptake mechanisms for all therapeutic ASOs.

3. **Endocytic uptake:** Receptor-mediated endocytosis via scavenger receptors (SR-A) is the dominant pathway, achieving 200246 nM in the phagolysosome. Macrophages have a 19x advantage over non-phagocytic cells.

4. **Clinical comparison:** Unlike mipomersen and inotersen, which require endosomal escape (~1-2% efficiency) to reach cytoplasmic mRNA, MRL-ASO-001's target (SL RNA) is in the phagolysosome — the SAME compartment where ASOs naturally accumulate after endocytosis. This eliminates the most significant bottleneck in ASO delivery.

**The PS backbone provides natural tropism for macrophages, making MRL-ASO-001 uniquely suited for anti-leishmanial therapy without complex formulation.**

## References

1. Geary RS et al. (2015) Adv Drug Deliv Rev 87:46-51 — PS-ASO pharmacokinetics
2. Crooke ST et al. (2017) Nat Rev Drug Discov 16(8):547-563 — ASO therapeutics review
3. Stein CA et al. (2010) Nucleic Acids Res 38(5):e3 — Gymnosis demonstration
4. Butler M et al. (2000) J Pharmacol Exp Ther 292(2):547-555 — SR-A binding of PS-ASOs
5. Parsegian VA (1969) Nature 221(5183):844-846 — Born solvation energy model
6. McLaughlin S (1989) Annu Rev Biophys Biophys Chem 18:113-136 — Membrane electrostatics
7. Steinman RM et al. (1976) J Cell Biol 68(3):665-687 — Macrophage pinocytosis
8. Koller E et al. (2011) Nucleic Acids Res 39(11):4795-4807 — ASO uptake in cells
9. Murphy MC et al. (2004) Biophys J 86(4):2530-2537 — ssDNA persistence length
10. Liang XH et al. (2015) Nucleic Acids Res 43(5):2927-2945 — PS-ASO cellular uptake
