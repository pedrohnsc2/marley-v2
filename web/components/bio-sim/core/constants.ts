import type { SceneMeta } from "./types";

export const BG_COLOR = 0x030a14;
export const FOG_DENSITY = 0.05;
export const STAR_COUNT = 800;
export const GROUP_ROTATION_SPEED = 0.0008;

export const STAGE_METAS: SceneMeta[] = [
  {
    name: "Subcutaneous Injection",
    shortName: "SC Injection",
    timeTag: "T = 0",
    color: "#93c5fd",
    description:
      "MRL-ASO-001 is administered via subcutaneous bolus injection. The 25-mer LNA-DNA-LNA gapmer with phosphorothioate backbone disperses into the extracellular matrix of the subcutaneous tissue.",
    metrics: [
      { label: "Bioavailability", value: "87.0%" },
      { label: "pH", value: "7.4" },
      { label: "Form", value: "Free ASO" },
    ],
  },
  {
    name: "Bloodstream Transit",
    shortName: "Bloodstream",
    timeTag: "T = 1–4 h",
    color: "#f87171",
    description:
      "ASO enters systemic circulation and binds plasma proteins (>85% bound). PS-backbone ensures nuclease resistance and extended half-life. Peak plasma concentration reached at ~3 hours.",
    metrics: [
      { label: "Cmax", value: "~2.4 μg/mL" },
      { label: "Tmax", value: "~3 h" },
      { label: "t½", value: "~21 days" },
    ],
  },
  {
    name: "Spleen & Liver Accumulation",
    shortName: "Organ Uptake",
    timeTag: "T = 4–12 h",
    color: "#fb923c",
    description:
      "PS-ASO accumulates preferentially in reticuloendothelial organs. Hepatic sinusoidal Kupffer cells and splenic macrophages capture circulating ASO via scavenger receptor-mediated uptake.",
    metrics: [
      { label: "Kp spleen", value: "~35×" },
      { label: "Kp liver", value: "~18×" },
      { label: "Uptake", value: "Scavenger-A" },
    ],
  },
  {
    name: "Macrophage Endocytosis",
    shortName: "Endocytosis",
    timeTag: "T = 12–24 h",
    color: "#c084fc",
    description:
      "Tissue-resident macrophages internalize ASO via clathrin-mediated endocytosis through Scavenger-A and CR3 receptors. Trimannose conjugation improves uptake 9.7× over naked ASO.",
    metrics: [
      { label: "Uptake", value: "CME / Scavenger" },
      { label: "EE pH", value: "6.5" },
      { label: "Vesicle", value: "~200 nm" },
    ],
  },
  {
    name: "Endosomal Trafficking",
    shortName: "Trafficking",
    timeTag: "T = 24–36 h",
    color: "#34d399",
    description:
      "ASO transits through the endosomal pathway: early endosome (pH 6.5) → multivesicular body (pH 5.5) → late endosome (pH 5.0). LNA modifications protect against lysosomal nucleases.",
    metrics: [
      { label: "pH gradient", value: "6.5 → 5.0 → 4.5" },
      { label: "Escape", value: "~1–2%" },
      { label: "Mechanism", value: "MVB retrofusion" },
    ],
  },
  {
    name: "Endosomal Escape & SL RNA Binding",
    shortName: "SL RNA Binding",
    timeTag: "T = 36–60 h",
    color: "#fbbf24",
    description:
      "1–2% of ASO escapes the endosome and reaches the cytoplasm. The escaped ASO hybridizes to the Spliced Leader RNA with high affinity (Kd ~2 nM), recruiting RNase H1 for catalytic cleavage.",
    metrics: [
      { label: "Kd", value: "~2 nM" },
      { label: "Target", value: "SL RNA (unique)" },
      { label: "Mechanism", value: "RNase H1" },
    ],
  },
  {
    name: "Leishmania Death",
    shortName: "Parasite Death",
    timeTag: "T = 60–72 h",
    color: "#e879f9",
    description:
      "RNase H1 cleavage of SL RNA abolishes trans-splicing in L. infantum, halting all mRNA maturation. The amastigote loses viability within the parasitophorous vacuole. Resistance is mathematically impossible.",
    metrics: [
      { label: "Reduction", value: ">90%" },
      { label: "EC50", value: "0.1 μM" },
      { label: "Resistance", value: "∞ years" },
    ],
  },
  {
    name: "TLR9 / Th1 Immune Activation",
    shortName: "TLR9 / Th1",
    timeTag: "T = 72 h+",
    color: "#4ade80",
    description:
      "CpG motifs in the PS-backbone activate endosomal TLR9, triggering NF-κB signaling and sustained Th1 cytokine release (IFN-γ, TNF-α). Dual mechanism: direct kill + immune priming.",
    metrics: [
      { label: "Pathway", value: "TLR9 → NF-κB" },
      { label: "Cytokines", value: "IFN-γ + TNF-α" },
      { label: "Response", value: "Th1 sustained" },
    ],
  },
];
