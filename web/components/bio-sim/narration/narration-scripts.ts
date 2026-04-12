import type { NarrationSegment } from "../core/types";

export const NARRATION_SCRIPTS: NarrationSegment[][] = [
  // Scene 0 — SC Injection
  [
    {
      text: "This is Marley — a computational journey to defeat visceral leishmaniasis. MRL-ASO-001 enters the body through subcutaneous injection, the same route used for insulin and many biologics.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "This twenty-five nucleotide antisense oligonucleotide carries an LNA-DNA-LNA gapmer design with a fully phosphorothioate backbone — engineered for nuclease resistance and tissue penetration.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "The drug disperses radially through the extracellular matrix, entering capillaries and lymphatic vessels. With eighty-seven percent bioavailability, nearly the entire dose reaches systemic circulation.",
      startPct: 0.60,
      durationPct: 0.28,
    },
  ],

  // Scene 1 — Bloodstream
  [
    {
      text: "In the bloodstream, MRL-ASO-001 rapidly binds plasma proteins — over eighty-five percent becomes protein-bound within minutes. This is not a limitation; it is the delivery mechanism.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "Protein binding shields the molecule from renal filtration, extending the terminal half-life to twenty-one days. The phosphorothioate backbone is the key: it creates high-affinity interactions with albumin.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "Peak plasma concentration reaches approximately 2.4 micrograms per milliliter at three hours post-injection. From here, the drug distributes preferentially to tissues rich in phagocytic cells.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],

  // Scene 2 — Organ Uptake
  [
    {
      text: "The liver and spleen become the primary reservoirs. Phosphorothioate ASOs accumulate in reticuloendothelial organs at concentrations thirty-five times the plasma level in the spleen.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "Kupffer cells lining the hepatic sinusoids and splenic macrophages capture circulating ASO through scavenger receptor class A — a receptor pathway that Leishmania-infected macrophages express at elevated levels.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "This is pharmacological serendipity: the very cells harboring L. infantum amastigotes are the cells that most avidly internalize the drug. The parasite's hiding place becomes its trap.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],

  // Scene 3 — Endocytosis
  [
    {
      text: "At the macrophage surface, scavenger receptors and complement receptor CR3 recognize the phosphorothioate backbone. Clathrin-coated pits form, initiating receptor-mediated endocytosis.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "Trimannose conjugation — a sugar tag recognized by mannose receptors on macrophages — improves cellular uptake by nearly ten-fold. The intracellular concentration reaches over two hundred thousand nanomolar.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "A pseudopod extends from the cell membrane, engulfing the ASO-receptor complex into an early endosomal vesicle approximately two hundred nanometers in diameter.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],

  // Scene 4 — Endosomal Trafficking
  [
    {
      text: "Inside the cell, the ASO enters the endosomal pathway — a gauntlet of progressively acidifying compartments. The early endosome maintains pH 6.5; the multivesicular body drops to 5.5.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "As pH falls, most therapeutic oligonucleotides are degraded. But the LNA modifications flanking the DNA gap protect MRL-ASO-001. The gapmer half-life exceeds one thousand hours — forty-five days of intracellular stability.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "The late endosome reaches pH 5.0. Here, a critical bottleneck: only one to two percent of internalized ASO molecules escape into the cytoplasm. The rest are routed to lysosomes.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],

  // Scene 5 — Endosomal Escape + SL RNA Binding (CLIMAX)
  [
    {
      text: "This is the decisive moment. Through multivesicular body back-fusion, a small fraction of ASO molecules breach the endosomal membrane and reach the cytoplasm — free to find their target.",
      startPct: 0.0,
      durationPct: 0.25,
    },
    {
      text: "The target is unique in all of biology: the Spliced Leader RNA. Every protein-coding transcript in Leishmania requires trans-splicing of this thirty-nine nucleotide leader sequence. There is no redundancy. No backup. No escape.",
      startPct: 0.27,
      durationPct: 0.25,
    },
    {
      text: "MRL-ASO-001 hybridizes to the SL RNA with a dissociation constant of two nanomolar and a binding energy of minus twenty-eight kilocalories per mole. The duplex recruits RNase H1, which catalytically cleaves the RNA strand.",
      startPct: 0.54,
      durationPct: 0.25,
    },
    {
      text: "One ASO molecule can destroy hundreds of SL RNA copies. Trans-splicing halts. mRNA maturation ceases. The parasite begins to die.",
      startPct: 0.82,
      durationPct: 0.16,
    },
  ],

  // Scene 6 — Parasite Death
  [
    {
      text: "Without functional mRNA, the amastigote cannot maintain its proteome. Membrane integrity fails. The cell shrinks within its parasitophorous vacuole as essential proteins are depleted.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "Computational modeling predicts greater than ninety percent parasite reduction. The EC50 is just 0.1 micromolar — a concentration easily achieved given the two hundred thousand nanomolar intracellular accumulation.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "And resistance? Mathematically impossible. Zero viable escape mutations exist that could disrupt ASO binding while preserving trans-splicing function. The target is evolutionarily locked.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],

  // Scene 7 — TLR9 / Th1 Immune Activation
  [
    {
      text: "But MRL-ASO-001 does more than kill parasites directly. The CpG dinucleotides in its phosphorothioate backbone activate Toll-like receptor 9 in the endosomal membrane.",
      startPct: 0.0,
      durationPct: 0.28,
    },
    {
      text: "TLR9 activation triggers NF-kappa-B signaling, unleashing a cascade of Th1 cytokines: interferon-gamma and TNF-alpha. This is precisely the immune profile needed to clear intracellular Leishmania.",
      startPct: 0.30,
      durationPct: 0.28,
    },
    {
      text: "The digital twin simulation confirms it: vaccine plus ASO dual therapy achieves complete parasite clearance in forty-one days. Vaccine alone never clears. This dual mechanism — direct kill plus immune priming — is the future of leishmaniasis therapy.",
      startPct: 0.62,
      durationPct: 0.28,
    },
  ],
];
