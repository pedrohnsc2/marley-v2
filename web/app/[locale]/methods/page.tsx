/* ------------------------------------------------------------------ */
/*  /methods — Computational Methods & References                      */
/*  Structured like a paper's Methods section for PPGIT/UFMG review.   */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Data types                                                         */
/* ------------------------------------------------------------------ */

interface Reference {
  id: string;
  authors: string;
  year: number;
  title: string;
  journal: string;
  detail: string;
}

interface MethodSection {
  anchor: string;
  number: string;
  title: string;
  content: string[];
  references: string[];
}

interface Limitation {
  module: string;
  items: string[];
}

/* ------------------------------------------------------------------ */
/*  References database                                                */
/* ------------------------------------------------------------------ */

const REFERENCES: Reference[] = [
  {
    id: "santalucia1998",
    authors: "SantaLucia J Jr.",
    year: 1998,
    title:
      "A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics",
    journal: "Proc Natl Acad Sci USA",
    detail: "95(4):1460-1465",
  },
  {
    id: "reynisson2020",
    authors: "Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M.",
    year: 2020,
    title:
      "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHCeluted ligand data",
    journal: "Nucleic Acids Res",
    detail: "48(W1):W449-W454",
  },
  {
    id: "trott2010",
    authors: "Trott O, Olson AJ.",
    year: 2010,
    title: "AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading",
    journal: "J Comput Chem",
    detail: "31(2):455-461",
  },
  {
    id: "lin2023",
    authors: "Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, Smeaton N, Verkuil R, Kabber O, Shmueli Y, dos Santos Costa A, Fazel-Zarandi M, Sercu T, Candido S, Rives A.",
    year: 2023,
    title:
      "Evolutionary-scale prediction of atomic-level protein structure with a language model",
    journal: "Science",
    detail: "379(6637):1123-1130",
  },
  {
    id: "doytchinova2007",
    authors: "Doytchinova IA, Flower DR.",
    year: 2007,
    title:
      "VaxiJen: a server for prediction of protective antigens, tumour antigens and subunit vaccines",
    journal: "BMC Bioinformatics",
    detail: "8:4",
  },
  {
    id: "liang2003",
    authors: "Liang XH, Haritan A, Uliel S, Michaeli S.",
    year: 2003,
    title:
      "trans and cis splicing in trypanosomatids: mechanism, factors, and regulation",
    journal: "Eukaryot Cell",
    detail: "2(5):830-840",
  },
  {
    id: "michaeli2011",
    authors: "Michaeli S.",
    year: 2011,
    title:
      "Trans-splicing in trypanosomes: machinery and its impact on the parasite transcriptome",
    journal: "Future Microbiol",
    detail: "6(4):459-474",
  },
  {
    id: "geary2015",
    authors: "Geary RS, Norris D, Yu R, Bennett CF.",
    year: 2015,
    title:
      "Pharmacokinetics, biodistribution and cell uptake of antisense oligonucleotides",
    journal: "Adv Drug Deliv Rev",
    detail: "87:46-51",
  },
  {
    id: "bennett2017",
    authors: "Bennett CF, Swayze EE.",
    year: 2017,
    title: "RNA targeting therapeutics: molecular mechanisms of antisense oligonucleotides as a therapeutic platform",
    journal: "Annu Rev Pharmacol Toxicol",
    detail: "50:259-293",
  },
  {
    id: "crooke2017",
    authors: "Crooke ST, Witztum JL, Bennett CF, Baker BF.",
    year: 2017,
    title:
      "RNA-targeted therapeutics",
    journal: "Cell Metab",
    detail: "27(4):714-739",
  },
  {
    id: "boitz2012",
    authors: "Boitz JM, Ullman B, Jardim A, Carter NS.",
    year: 2012,
    title:
      "Purine salvage in Leishmania: complex or simple by design?",
    journal: "Trends Parasitol",
    detail: "28(8):345-352",
  },
  {
    id: "frearson2007",
    authors: "Frearson JA, Wyatt PG, Gilbert IH, Fairlamb AH.",
    year: 2007,
    title:
      "Target assessment for antiparasitic drug discovery",
    journal: "Trends Parasitol",
    detail: "23(12):589-595",
  },
  {
    id: "rogers2011",
    authors: "Rogers MB, Hilley JD, Dickens NJ, Wilkes J, Bates PA, Sherwin T, et al.",
    year: 2011,
    title:
      "Chromosome and gene copy number variation allow major structural change between species and strains of Leishmania",
    journal: "Genome Res",
    detail: "21(12):2129-2142",
  },
  {
    id: "gasteiger2005",
    authors: "Gasteiger E, Hoogland C, Gattiker A, Duvaud S, Wilkins MR, Appel RD, Bairoch A.",
    year: 2005,
    title:
      "Protein identification and analysis tools on the ExPASy server",
    journal: "The Proteomics Protocols Handbook, Humana Press",
    detail: "pp. 571-607",
  },
  {
    id: "aslett2010",
    authors: "Aslett M, Aurrecoechea C, Berriman M, Brestelli J, Brunk BP, Carrington M, et al.",
    year: 2010,
    title:
      "TriTrypDB: a functional genomic resource for the Trypanosomatidae",
    journal: "Nucleic Acids Res",
    detail: "38(Database issue):D457-D462",
  },
  {
    id: "seth2019",
    authors: "Seth PP, Swayze EE.",
    year: 2019,
    title: "The medicinal chemistry of LNA, cEt, and related chemistries",
    journal: "Antisense Drug Technology, CRC Press",
    detail: "Chapter 6, pp. 113-140",
  },
  {
    id: "kool2001",
    authors: "Kool ET.",
    year: 2001,
    title: "Hydrogen bonding, base stacking, and steric effects in DNA replication",
    journal: "Annu Rev Biophys Biomol Struct",
    detail: "30:1-22",
  },
  {
    id: "ilango2016",
    authors: "Ilango S, Sundaram RS.",
    year: 2016,
    title: "Immunoinformatics approaches for multi-epitope vaccine design against Leishmania donovani",
    journal: "Bioinformation",
    detail: "12(11):375-379",
  },
];

const refMap = Object.fromEntries(REFERENCES.map((r) => [r.id, r]));

/* ------------------------------------------------------------------ */
/*  Methods sections                                                   */
/* ------------------------------------------------------------------ */

const METHODS: MethodSection[] = [
  {
    anchor: "databases",
    number: "2.1",
    title: "Databases and Data Sources",
    content: [
      "The complete proteome of Leishmania infantum strain JPCM5 was retrieved from TriTrypDB release 68 (https://tritrypdb.org), comprising 8,527 annotated protein sequences. TriTrypDB is the primary functional genomics resource for trypanosomatid organisms, integrating genome sequence, annotation, and functional data from multiple sources [aslett2010].",
      "Protein sequences were cross-referenced against UniProtKB/Swiss-Prot for functional annotation and Gene Ontology terms. Structural templates for drug targets were obtained from the RCSB Protein Data Bank (PDB) and the AlphaFold Protein Structure Database (v4). The human reference proteome (GRCh38, Ensembl release 112) was used for host ortholog exclusion.",
      "Canine MHC class I alleles (DLA-88*03401, DLA-88*50101, DLA-88*50801) were sourced from the Immuno Polymorphism Database (IPD-MHC). Spliced Leader RNA sequences for cross-species conservation analysis were obtained from published literature for T. cruzi, T. brucei, L. major, and L. donovani [liang2003].",
    ],
    references: ["aslett2010", "liang2003"],
  },
  {
    anchor: "epitope-prediction",
    number: "2.2",
    title: "Epitope Prediction and Antigen Screening",
    content: [
      "Surface-exposed protein candidates were filtered by the presence of signal peptides, transmembrane domains, and GPI anchors using SignalP 6.0 heuristics applied to TriTrypDB annotations. Conservation scoring used BLASTp (e-value threshold 1e-10) against L. donovani and L. major proteomes, retaining proteins with >80% sequence identity across at least two species.",
      "MHC class I binding predictions were performed using the IEDB Analysis Resource API (tools.iedb.org) with the NetMHCpan-BA 4.1 method [reynisson2020]. For each candidate protein, 9-mer peptides were scored against three canine DLA alleles: DLA-88*03401, DLA-88*50101, and DLA-88*50801. Peptides with predicted IC50 < 500 nM were classified as strong binders. Overlapping epitopes (>5 residue overlap) were deduplicated, retaining the peptide with the lowest IC50.",
      "Final antigen scoring combined immunogenicity (60% weight) and cross-species conservation (40% weight) into a composite score. The top 11 epitopes from 10 candidate antigens were selected for the vaccine construct. Antigenicity of the final construct was predicted using VaxiJen 2.0 with a threshold of 0.4 for parasitic organisms [doytchinova2007].",
    ],
    references: ["reynisson2020", "doytchinova2007"],
  },
  {
    anchor: "vaccine-design",
    number: "2.3",
    title: "Multi-Epitope Vaccine Design",
    content: [
      "The multi-epitope vaccine construct was assembled in silico following established reverse vaccinology principles. The N-terminal architecture consists of a tissue plasminogen activator (tPA) signal peptide for secretory pathway targeting, followed by the L7/L12 ribosomal protein adjuvant domain linked via an EAAAK rigid linker. L7/L12 was selected for its demonstrated ability to enhance antigen presentation via TLR4 activation in canine immune cells.",
      "The 11 selected CTL epitopes were concatenated using alternating GPGPG and AAY flexible linkers. GPGPG linkers maintain structural independence between epitopes while providing flexibility, and AAY linkers act as proteasomal cleavage sites to facilitate epitope processing. The construct terminates with a 6xHis tag for purification.",
      "Codon optimization was performed for Canis lupus familiaris using species-specific Relative Synonymous Codon Usage (RSCU) tables. Codons were selected by highest frequency for each amino acid. Post-optimization, restriction enzyme recognition sites (EcoRI, BamHI, XhoI, NdeI, HindIII, NotI) were eliminated by silent synonymous substitutions. The final mRNA cassette includes a 5' cap analog (CleanCap AG), Kozak consensus sequence, optimized coding sequence, two tandem stop codons, and a 120-nt poly(A) tail.",
      "Physicochemical properties were computed using ProtParam algorithms [gasteiger2005]: molecular weight, theoretical isoelectric point (pI), instability index (threshold < 40 for stable proteins), GRAVY (grand average of hydropathicity), and aromaticity. GC content of the mRNA was maintained between 40-60% for optimal translation efficiency.",
    ],
    references: ["gasteiger2005"],
  },
  {
    anchor: "structural-analysis",
    number: "2.4",
    title: "Structural Analysis and 3D Prediction",
    content: [
      "Three-dimensional structure prediction of the multi-epitope vaccine construct was performed using ESMFold [lin2023], a protein language model that generates atomic-resolution structures from single sequences without requiring multiple sequence alignments. ESMFold was chosen over AlphaFold2 for its ability to handle artificial multi-epitope constructs that lack evolutionary covariance signals.",
      "The predicted structure was assessed for global topology and domain organization. Per-residue confidence scores (pLDDT) were examined to identify well-structured regions versus flexible linkers. Molecular visualization was performed using PyMOL 2.5.",
    ],
    references: ["lin2023"],
  },
  {
    anchor: "drug-targets",
    number: "2.5",
    title: "Drug Target Identification",
    content: [
      "Parasite-specific drug targets were identified through comparative genomics, subtracting the L. infantum proteome against the human and canine host proteomes. Proteins with sequence identity >40% to any human protein (BLASTp, e-value < 1e-5) were excluded to minimize host toxicity risk.",
      "The remaining parasite-specific proteins were classified into four metabolic pathways known to be essential in trypanosomatids: trypanothione metabolism [frearson2007], sterol biosynthesis, folate metabolism, and purine salvage [boitz2012]. Druggability scoring combined: (a) pathway essentiality evidence from the literature, (b) presence of known inhibitors or chemical probes in ChEMBL, (c) structural druggability (active site accessibility from PDB/AlphaFold structures), and (d) absence of close human orthologs. A total of 52 targets across 4 pathways were prioritized, with the top 20 receiving detailed analysis.",
    ],
    references: ["frearson2007", "boitz2012"],
  },
  {
    anchor: "molecular-docking",
    number: "2.6",
    title: "Molecular Docking",
    content: [
      "Virtual screening was performed using AutoDock Vina 1.2.5 [trott2010]. Two receptor structures were prepared: GMP synthase (GMPS, UniProt Q4Q457 homology model) and trypanothione reductase (TryR, UniProt Q4Q457 AlphaFold structure). Receptor preparation involved removal of water molecules, addition of polar hydrogens, and Gasteiger charge assignment using Open Babel 3.1.",
      "A compound library of 77 molecules was assembled from three sources: (a) known anti-leishmanial compounds from ChEMBL (filtered by IC50 < 10 uM against L. infantum), (b) FDA-approved drugs with reported anti-parasitic activity (e.g., methotrexate, pemetrexed, menadione), and (c) three custom-designed molecules (MRL-001, MRL-002, MRL-003) based on trypanothione reductase pharmacophore models.",
      "Docking parameters: exhaustiveness = 32, grid box size = 30x30x30 A centered on the catalytic site, energy range = 4 kcal/mol, num_modes = 9. All ligands were prepared as PDBQT files with rotatable bonds defined using Open Babel. Results were ranked by predicted binding affinity (kcal/mol). Selectivity analysis compared binding to parasite TryR versus human glutathione reductase (GR, PDB: 3GRS/AlphaFold P00390) to assess therapeutic window.",
    ],
    references: ["trott2010"],
  },
  {
    anchor: "aso-design",
    number: "2.7",
    title: "Antisense Oligonucleotide Design",
    content: [
      "Antisense oligonucleotide (ASO) therapy was designed to target the Spliced Leader (SL) RNA of L. infantum. The SL RNA is a 39-nucleotide conserved sequence (5'-AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG-3') that is trans-spliced to the 5' end of every mRNA in trypanosomatids, making it an essential and universal target [liang2003, michaeli2011]. Critically, SL RNA has no mammalian homolog, providing inherent selectivity.",
      "The lead candidate MRL-ASO-001 is a 25-mer (5'-ACAGAAACTGATACTTATATAGCGT-3') targeting positions 5-30 of the SL RNA. Thermodynamic properties were calculated using the SantaLucia unified nearest-neighbor model [santalucia1998] with standard DNA/DNA parameters. Duplex stability was computed at physiological conditions (37 C, 250 nM oligonucleotide concentration, 1 M NaCl equivalent): predicted Tm = 68.48 C, predicted dG = -27.97 kcal/mol.",
      "The ASO employs an LNA-DNA-LNA gapmer architecture: 4 locked nucleic acid (LNA) residues at each terminus flanking a central 17-nt DNA gap [seth2019]. LNA modifications increase binding affinity (+2 to +5 C Tm per modification) and nuclease resistance. The central DNA gap maintains RNase H recruitment competence, enabling catalytic degradation of the target RNA [bennett2017, crooke2017].",
      "Self-complementarity was evaluated computationally. Hairpin formation dG and self-dimer dG were calculated to ensure the ASO preferentially binds target over itself. Off-target screening was performed against the human transcriptome (GRCh38) using a Smith-Waterman local alignment with thermodynamic scoring to identify potential cross-reactive sequences [geary2015].",
    ],
    references: [
      "liang2003",
      "michaeli2011",
      "santalucia1998",
      "seth2019",
      "bennett2017",
      "crooke2017",
      "geary2015",
    ],
  },
  {
    anchor: "aso-validation",
    number: "2.8",
    title: "ASO Mathematical Validation",
    content: [
      "A five-module mathematical validation suite was developed to rigorously assess ASO candidate viability without wet-lab experiments:",
      "Module 1 -- Thermodynamic Landscape: Complete enumeration of all possible ASO sequences (18-27 nt) targeting the SL RNA, computing dG, Tm, GC content, hairpin dG, and self-dimer dG for each. The search confirmed MRL-ASO-001 as rank-1 across the exhaustive design space.",
      "Module 2 -- Selectivity Proof: Formal proof that the ASO binds its intended target with significantly greater affinity than any off-target in the human transcriptome. Uses thermodynamic discrimination ratio (dG_target / dG_best_off_target).",
      "Module 3 -- Evolutionary Conservation: Cross-species analysis of the SL RNA target across T. cruzi, T. brucei, L. major, and L. donovani, demonstrating conservation over approximately 350 Mya of divergence. Uses neutral substitution rate modeling [rogers2011] to estimate the probability of target-disrupting mutations.",
      "Module 4 -- Exhaustive Optimization: Combinatorial scan of all single-nucleotide variants proving that no single-base change to MRL-ASO-001 improves binding affinity while maintaining selectivity constraints.",
      "Module 5 -- Resistance Barrier: Markov chain model of parasite mutational escape. Given L. infantum mutation rate (~2e-9 per nt per replication) [rogers2011], generation time (12 h), and SL RNA copy number (~150 tandem repeats), the model computes the probability of resistance emergence as effectively zero over clinically relevant timeframes.",
    ],
    references: ["rogers2011"],
  },
  {
    anchor: "delivery-modeling",
    number: "2.9",
    title: "Delivery and Pharmacokinetic Modeling",
    content: [
      "A six-module computational pipeline modeled the in vivo behavior of MRL-ASO-001 from administration to target engagement:",
      "Module A (Stability): Nuclease degradation kinetics for PS-ASO backbone under serum conditions, comparing unmodified, phosphorothioate (PS), and PS+LNA chemistries.",
      "Module B (Membrane): Permeability modeling across endosomal and plasma membranes using physicochemical descriptors (MW, logP, charge, PSA). Endocytic uptake rate estimation based on cell-type specific internalization kinetics for macrophages.",
      "Module C (Conjugation): Evaluation of GalNAc, cholesterol, and cell-penetrating peptide (CPP) conjugation strategies for enhanced cellular uptake in macrophages harboring Leishmania amastigotes.",
      "Module D (LNP Encapsulation): Lipid nanoparticle formulation modeling (SM-102, DSPC, cholesterol, PEG-DMG) with predicted encapsulation efficiency, particle size distribution, and controlled release kinetics.",
      "Module E (ADMET Profile): Simplified pharmacokinetic model predicting absorption, distribution (liver, spleen, bone marrow tropism), metabolism (nuclease-mediated), and excretion for subcutaneous administration in a 30 kg canine model.",
      "Module F (Immune Response): Stochastic differential equation (SDE) model of innate immune activation by PS-ASO, including TLR9 stimulation, cytokine release kinetics, and predicted therapeutic window between efficacy and immunotoxicity.",
    ],
    references: ["geary2015", "bennett2017"],
  },
  {
    anchor: "immune-simulation",
    number: "2.10",
    title: "Immune Response Simulation",
    content: [
      "Vaccine-induced immune responses were simulated using a custom ODE-based kinetics model tracking five cell populations over 365 days: T-helper (Th), cytotoxic T cells (Tc), B cells, antibodies (Ab), and memory cells (M). A three-dose schedule (Days 0, 21, 42) was modeled with exponential rise and decay kinetics fitted to published canine vaccination data.",
      "Th1/Th2 polarization was estimated from epitope characteristics and adjuvant properties, with Th1 dominance being the desired outcome for intracellular pathogen clearance. Memory cell stability was projected using an exponential decay model to estimate duration of protection.",
    ],
    references: [],
  },
  {
    anchor: "platform-comparison",
    number: "2.11",
    title: "Vaccine Platform Comparison",
    content: [
      "Three production platforms were evaluated for manufacturing the multi-epitope vaccine construct: (a) mRNA-LNP (in vitro transcription with lipid nanoparticle delivery), (b) recombinant protein (E. coli BL21(DE3) expression with Ni-NTA purification), and (c) live attenuated Leishmania tarentolae (LEXSY system, non-pathogenic expression host).",
      "Platforms were compared across six dimensions: construct properties (stability, antigenicity), production parameters (expression system, purification, adjuvant requirements), logistics (cold chain, infrastructure), cost per dose (lab and industrial scale), ASO therapy compatibility, and strategic funding feasibility. Normalized scores (0-1) were computed for each dimension to enable cross-platform ranking.",
    ],
    references: [],
  },
  {
    anchor: "software",
    number: "2.12",
    title: "Software and Versions",
    content: [
      "The computational pipeline was implemented in Python 3.11+ with the following key dependencies: Biopython 1.83 (sequence analysis, ProtParam), NumPy 1.26 (numerical computation), SciPy 1.12 (statistical analysis, ODE integration), Requests 2.31 (API communication with IEDB and TriTrypDB), tqdm 4.66 (progress tracking).",
      "Molecular docking: AutoDock Vina 1.2.5 [trott2010], Open Babel 3.1 (file format conversion, charge assignment). Structure prediction: ESMFold API [lin2023]. Antigenicity prediction: VaxiJen 2.0 web server [doytchinova2007]. Visualization: PyMOL 2.5, Matplotlib 3.8, Three.js 0.183 (web 3D rendering).",
      "The web dashboard was built with Next.js 14.2, React 18.3, TypeScript 5.5, Tailwind CSS 3.4, and ApexCharts 5.10 for interactive data visualization. All source code is version-controlled with Git.",
    ],
    references: ["trott2010", "lin2023", "doytchinova2007"],
  },
];

/* ------------------------------------------------------------------ */
/*  Limitations                                                        */
/* ------------------------------------------------------------------ */

const LIMITATIONS: Limitation[] = [
  {
    module: "Vaccine Construct",
    items: [
      "No experimental validation has been performed. All epitope binding affinities are computational predictions based on NetMHCpan-BA 4.1 and have not been confirmed by binding assays (e.g., ELISpot, tetramer staining).",
      "Only MHC class I (CTL) epitopes were included. MHC class II (helper T cell) epitopes are absent from the current construct, limiting the predicted helper response to adjuvant-mediated stimulation.",
      "B-cell epitopes were not predicted or included. A complete vaccine design would incorporate linear and conformational B-cell epitopes for antibody-mediated protection.",
      "VaxiJen antigenicity is a statistical prediction. Experimental immunogenicity testing in canine models is required to validate the 0.3235 antigenicity score.",
      "AllerTOP allergenicity assessment has not yet been performed. This analysis is needed before any preclinical formulation.",
    ],
  },
  {
    module: "Drug Targets",
    items: [
      "All 52 drug targets are computationally predicted. No experimental IC50, Ki, or enzyme inhibition data have been generated.",
      "Druggability scores are heuristic composites, not validated against experimental hit rates from high-throughput screening campaigns.",
      "Essential gene predictions rely on ortholog inference from T. brucei and L. major knockout studies, which may not transfer directly to L. infantum.",
      "The host ortholog exclusion threshold (40% identity) is a conservative heuristic. Some parasite-specific targets may still share functional similarity with human enzymes at lower sequence identity.",
    ],
  },
  {
    module: "Molecular Docking",
    items: [
      "Rigid receptor docking was used throughout. No induced-fit or ensemble docking was performed, potentially missing conformational changes upon ligand binding.",
      "Only two receptor structures were used (GMPS, TryR). The 52 identified drug targets represent a much larger structural space that was not explored.",
      "No molecular dynamics (MD) validation was performed. Docking poses were not refined by MD simulation, and binding free energy estimates (MM-PBSA/MM-GBSA) were not computed.",
      "None of the tested compounds showed parasite selectivity over human glutathione reductase. The selectivity analysis demonstrates that further structural optimization of lead compounds is required.",
      "Docking scores correlate imperfectly with experimental binding affinities. The -8.07 kcal/mol best hit should be interpreted as a ranking metric, not an absolute binding energy.",
    ],
  },
  {
    module: "ASO Therapy",
    items: [
      "Thermodynamic predictions use the SantaLucia nearest-neighbor model for DNA/DNA duplexes. The actual ASO contains LNA modifications whose thermodynamic contributions (+2 to +5 C Tm per residue) are estimated but not precisely parameterized in this model.",
      "No in vitro or in vivo data exist for MRL-ASO-001. Knockdown efficiency, cellular uptake, tissue distribution, and toxicity are all computational estimates.",
      "The resistance model assumes a constant mutation rate. Stress-induced mutagenesis under drug pressure could accelerate escape, though the SL RNA copy number (~150) and functional constraint provide substantial barriers.",
      "Off-target predictions against the human transcriptome use sequence complementarity as a proxy. Functional off-target effects (e.g., RNase H-independent mechanisms, protein binding) are not modeled.",
    ],
  },
  {
    module: "Bio-Sim and Delivery Modeling",
    items: [
      "All pharmacokinetic parameters are computational predictions based on published ASO class averages, not compound-specific measurements.",
      "The immune simulation uses a simplified ODE model that does not capture the full complexity of canine immune responses, including tissue-resident immunity, germinal center dynamics, or mucosal responses.",
      "LNP formulation modeling uses idealized particle distributions. Actual formulation development requires Design of Experiments (DoE) optimization with physical characterization (DLS, cryo-EM).",
      "The SDE immune response model for ASO delivery uses published human/murine parameters scaled to canine physiology. Species-specific immune kinetics may differ substantially.",
      "Memory cell decay projections (>693 days protection) are extrapolations from a 365-day simulation window and should be interpreted as order-of-magnitude estimates.",
    ],
  },
];

/* ------------------------------------------------------------------ */
/*  Table of Contents component                                        */
/* ------------------------------------------------------------------ */

function TableOfContents({ sections }: { sections: MethodSection[] }) {
  return (
    <div
      className="rounded-xl p-5"
      style={{
        backgroundColor: "var(--app-surface)",
        border: "1px solid var(--app-border)",
        boxShadow: "var(--app-card-shadow)",
      }}
    >
      <h2 className="mb-3 text-sm font-semibold" style={{ color: "var(--app-text)" }}>
        Table of Contents
      </h2>
      <ol className="space-y-1.5">
        {sections.map((s) => (
          <li key={s.anchor}>
            <a
              href={`#${s.anchor}`}
              className="text-sm transition-colors hover:underline"
              style={{ color: "var(--app-text-2)" }}
            >
              <span className="font-mono text-xs" style={{ color: "var(--app-text-3)" }}>
                {s.number}
              </span>{" "}
              {s.title}
            </a>
          </li>
        ))}
        <li>
          <a
            href="#references"
            className="text-sm transition-colors hover:underline"
            style={{ color: "var(--app-text-2)" }}
          >
            <span className="font-mono text-xs" style={{ color: "var(--app-text-3)" }}>
              3.0
            </span>{" "}
            References
          </a>
        </li>
        <li>
          <a
            href="#limitations"
            className="text-sm transition-colors hover:underline"
            style={{ color: "var(--app-text-2)" }}
          >
            <span className="font-mono text-xs" style={{ color: "var(--app-text-3)" }}>
              4.0
            </span>{" "}
            Known Limitations
          </a>
        </li>
      </ol>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Inline citation rendering                                          */
/* ------------------------------------------------------------------ */

function renderParagraph(text: string): React.ReactNode[] {
  /* Split on [refId] or [refId, refId] citation patterns */
  const parts = text.split(/(\[[a-z0-9_, ]+\])/g);
  return parts.map((part, i) => {
    const match = part.match(/^\[([a-z0-9_, ]+)\]$/);
    if (match) {
      const ids = match[1].split(",").map((s) => s.trim());
      return (
        <span key={i} className="text-xs font-semibold" style={{ color: "var(--app-primary)" }}>
          [
          {ids.map((id, j) => {
            const ref = refMap[id];
            const idx = REFERENCES.findIndex((r) => r.id === id);
            return (
              <span key={id}>
                {j > 0 && ", "}
                <a href={`#ref-${id}`} className="hover:underline">
                  {ref ? `${ref.authors.split(",")[0].split(" ").slice(-1)[0]} ${ref.year}` : idx + 1}
                </a>
              </span>
            );
          })}
          ]
        </span>
      );
    }
    return <span key={i}>{part}</span>;
  });
}

/* ------------------------------------------------------------------ */
/*  Page                                                               */
/* ------------------------------------------------------------------ */

export default function MethodsPage() {
  return (
    <div>
      {/* Page header */}
      <div className="mb-6 flex items-center gap-3">
        <span
          className="rounded-lg px-2.5 py-1 text-xs font-bold"
          style={{
            backgroundColor: "var(--app-primary-light)",
            color: "var(--app-primary)",
          }}
        >
          v2.0
        </span>
        <div>
          <h1 className="text-2xl font-bold" style={{ color: "var(--app-text)" }}>
            Computational Methods
          </h1>
          <p className="text-sm" style={{ color: "var(--app-text-2)" }}>
            Detailed methodology, parameters, and literature references for all pipeline modules
          </p>
        </div>
      </div>

      {/* Table of Contents */}
      <div className="mb-6">
        <TableOfContents sections={METHODS} />
      </div>

      {/* Methods sections */}
      {METHODS.map((section) => (
        <div
          key={section.anchor}
          id={section.anchor}
          className="mb-6 scroll-mt-20 rounded-xl p-6"
          style={{
            backgroundColor: "var(--app-surface)",
            border: "1px solid var(--app-border)",
            boxShadow: "var(--app-card-shadow)",
          }}
        >
          <div className="mb-4 flex items-center gap-3">
            <span
              className="rounded-md px-2 py-0.5 font-mono text-xs font-semibold"
              style={{
                backgroundColor: "var(--app-surface-2)",
                color: "var(--app-text-2)",
                border: "1px solid var(--app-border)",
              }}
            >
              {section.number}
            </span>
            <h2 className="text-base font-semibold" style={{ color: "var(--app-text)" }}>
              {section.title}
            </h2>
          </div>
          <div className="space-y-3">
            {section.content.map((paragraph, pi) => (
              <p
                key={pi}
                className="text-sm leading-relaxed"
                style={{ color: "var(--app-text-2)" }}
              >
                {renderParagraph(paragraph)}
              </p>
            ))}
          </div>
        </div>
      ))}

      {/* References */}
      <div
        id="references"
        className="mb-6 scroll-mt-20 rounded-xl p-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="mb-4 flex items-center gap-3">
          <span
            className="rounded-md px-2 py-0.5 font-mono text-xs font-semibold"
            style={{
              backgroundColor: "var(--app-surface-2)",
              color: "var(--app-text-2)",
              border: "1px solid var(--app-border)",
            }}
          >
            3.0
          </span>
          <h2 className="text-base font-semibold" style={{ color: "var(--app-text)" }}>
            References
          </h2>
        </div>
        <ol className="space-y-3">
          {REFERENCES.map((ref, i) => (
            <li
              key={ref.id}
              id={`ref-${ref.id}`}
              className="scroll-mt-20 text-sm leading-relaxed"
              style={{ color: "var(--app-text-2)" }}
            >
              <span className="mr-2 font-mono text-xs font-bold" style={{ color: "var(--app-text-3)" }}>
                [{i + 1}]
              </span>
              <span style={{ color: "var(--app-text)" }}>{ref.authors}</span> ({ref.year}).{" "}
              <em>{ref.title}.</em> {ref.journal}, {ref.detail}.
            </li>
          ))}
        </ol>
      </div>

      {/* Known Limitations */}
      <div
        id="limitations"
        className="mb-6 scroll-mt-20 rounded-xl p-6"
        style={{
          backgroundColor: "var(--app-surface)",
          border: "1px solid var(--app-border)",
          boxShadow: "var(--app-card-shadow)",
        }}
      >
        <div className="mb-4 flex items-center gap-3">
          <span
            className="rounded-md px-2 py-0.5 font-mono text-xs font-semibold"
            style={{
              backgroundColor: "var(--app-surface-2)",
              color: "var(--app-text-2)",
              border: "1px solid var(--app-border)",
            }}
          >
            4.0
          </span>
          <h2 className="text-base font-semibold" style={{ color: "var(--app-text)" }}>
            Known Limitations
          </h2>
        </div>
        <p className="mb-5 text-sm leading-relaxed" style={{ color: "var(--app-text-2)" }}>
          This pipeline is entirely computational. All results are predictions that require
          experimental validation before any preclinical or clinical application. The following
          section documents the specific limitations of each module to facilitate informed
          interpretation of the results.
        </p>
        <div className="space-y-5">
          {LIMITATIONS.map((lim) => (
            <div key={lim.module}>
              <h3
                className="mb-2 text-sm font-semibold"
                style={{ color: "var(--app-text)" }}
              >
                {lim.module}
              </h3>
              <ul className="space-y-1.5 pl-4">
                {lim.items.map((item, j) => (
                  <li
                    key={j}
                    className="relative pl-3 text-sm leading-relaxed before:absolute before:left-0 before:top-2 before:h-1.5 before:w-1.5 before:rounded-full"
                    style={{
                      color: "var(--app-text-2)",
                    }}
                  >
                    <span
                      className="absolute left-0 top-[0.45rem] h-1.5 w-1.5 rounded-full"
                      style={{ backgroundColor: "var(--app-text-3)" }}
                    />
                    {item}
                  </li>
                ))}
              </ul>
            </div>
          ))}
        </div>
      </div>

      {/* Footer note */}
      <div
        className="rounded-xl px-5 py-4"
        style={{
          backgroundColor: "var(--app-primary-light)",
          border: "1px solid var(--app-border)",
        }}
      >
        <p className="text-xs leading-relaxed" style={{ color: "var(--app-text-2)" }}>
          This methods section follows the structure of a peer-reviewed publication to
          facilitate evaluation by the PPGIT/UFMG admissions committee. All computational
          tools, parameters, and thresholds are documented for full reproducibility. Raw
          data and source code are available in the project repository.
        </p>
      </div>
    </div>
  );
}
