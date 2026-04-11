# How I Built a Drug Discovery Platform After Losing My Dog to Leishmaniasis

## 84,000 lines of code, 11 AI modules, and a molecule that might save other dogs' lives

---

In 2024, I lost my dog Marley to canine visceral leishmaniasis — a parasitic disease transmitted by sandflies that kills thousands of dogs every year in Brazil. There's no cure. The only approved vaccine (Leish-Tec) was suspended in 2023. Treatment options are toxic, expensive, and often ineffective.

I'm a software developer, not a biologist. But I know how to build systems, write code, and make computers solve problems. So I did the only thing I knew how to do: I started coding.

What began as grief turned into a computational drug discovery platform with 84,000+ lines of Python, 455 automated tests, and a therapeutic molecule candidate that — at least computationally — looks very promising.

This is the story of how that happened.

---

## The Starting Point: 8,527 Proteins

The *Leishmania infantum* parasite (the species that killed Marley) has a fully sequenced genome available in public databases. I downloaded all 8,527 annotated protein sequences and asked a simple question: **which of these proteins could be a good vaccine target?**

The filtering pipeline was straightforward:

- **Surface exposure**: only proteins the immune system can actually "see" (filtered with SignalP 6.0)
- **Conservation**: proteins that are the same across different Leishmania strains (if it mutates easily, it's a bad target)
- **Immunogenicity**: proteins that canine immune cells can recognize (predicted with IEDB tools against 3 dog MHC alleles)

8,527 proteins became 139 surface proteins, which became **11 epitopes** — short peptide fragments that the dog's immune system can recognize and attack, with binding affinities between 11 and 118 nanomolar.

These 11 epitopes became the foundation for everything that followed.

![Marley vaccine 3D structure](marley_vaccine_3d.png)
*The 3D structure of our multi-epitope vaccine construct (335 amino acids), predicted by ESMFold. Orange: signal peptide. Blue: adjuvant. Green: the 11 epitopes from L. infantum.*

---

## The Pivot: When Small Molecules Failed

With the vaccine designed, I moved to the second front: finding a drug that could treat already-infected dogs.

I mapped 52 enzymes across 5 metabolic pathways of the parasite, identified 7 priority drug targets (enzymes that exist in the parasite but not in humans), and ran molecular docking simulations with AutoDock Vina — essentially testing whether known drugs could fit into these enzyme "locks."

The results were mixed. Pentamidine (a known antileishmanial drug) ranked #2, which validated the methodology. But when I designed custom molecules (the MRL series), I hit a wall.

**MRL-003**, my best custom molecule, bound *better* to a human enzyme than to the parasite target (-8.68 vs -7.74 kcal/mol). I tried 20 redesigns. Zero were selective.

This was an honest negative result, and I documented it as such. But it forced a critical pivot: **if small molecules won't work, what about targeting RNA instead?**

![Docking visualization](gmps_methotrexate_docking.png)
*Molecular docking of methotrexate against the GMPS enzyme of L. infantum. The docking methodology was validated, but custom molecules failed the selectivity test.*

---

## Information Theory Meets Parasitology

Here's where things got interesting. I applied **Shannon entropy** — a concept from information theory, not biology — to the entire transcriptome of the parasite.

Shannon entropy measures how "predictable" something is. In DNA/RNA, a position where you always find the same nucleotide has low entropy (very predictable, very conserved). A position that varies randomly has high entropy (unpredictable, not important).

I computed entropy at every position across all the parasite's RNA, then asked: **what sequence is the most conserved in the parasite AND the most absent in humans?**

The answer was unambiguous: the **Spliced Leader RNA** (SL RNA).

The SL RNA is a 39-nucleotide sequence that gets attached to the beginning of *every single mRNA* the parasite makes. It's been conserved for approximately 500 million years of evolution. It's essential — without it, the parasite can't make any proteins. And here's the critical part: **it doesn't exist in humans or dogs.**

This is the dream target for a drug: absolutely essential for the parasite, completely absent from the host.

---

## Designing the Molecule: MRL-ASO-001

With the target identified, I designed an **antisense oligonucleotide (ASO)** — a short synthetic strand of modified DNA that binds to the SL RNA and blocks it.

Out of 119 candidate ASO sequences I evaluated, one stood out:

**MRL-ASO-001**: a 25-nucleotide gapmer with LNA (Locked Nucleic Acid) modifications on the ends and a phosphorothioate backbone throughout. It binds the SL RNA with a predicted energy of -27.97 kcal/mol (very strong), has zero off-targets in the human and dog genomes, and — because of its chemical design — has a dual mechanism:

1. **Antisense**: physically blocks the SL RNA, preventing all protein production in the parasite
2. **Immunostimulatory**: the phosphorothioate backbone activates TLR9 receptors on immune cells, triggering an innate immune response against the parasite

One molecule, two mechanisms. Block the parasite's biology AND wake up the dog's immune system.

![ASO 3D structure](mrl_aso_001_3d.png)
*3D structure of MRL-ASO-001. The molecule is a 25-nucleotide gapmer designed to bind the Spliced Leader RNA of L. infantum.*

---

## Mathematical Validation: Is This Molecule Actually Good?

Designing a molecule is one thing. Proving it's good is another. I built a mathematical validation suite with 6 independent dimensions, each testing a different aspect of the molecule.

### The Math Certificate: 52/60

**1. Thermodynamic Optimality (6/10)**
The binding energy is -27.97 kcal/mol. Good, but not the theoretical optimum — a 30-nucleotide version would bind even stronger. I deliberately chose 25 nt because shorter molecules are easier to deliver into cells. This is an honest trade-off.

**2. Information Geometry (10/10)**
Using Fisher-Rao distances (a concept from differential geometry), I measured how "mathematically alien" the SL RNA is compared to human and dog RNA. The result: the SL RNA sits **14.5 times farther** from host sequences than host sequences are from each other. In the geometry of information, this target is on another planet.

**3. Topological Stability (10/10)**
Using Topological Data Analysis (persistent homology), I showed that when MRL-ASO-001 binds the SL RNA, it creates 5 stable topological features in the resulting complex. The binding doesn't just happen — it creates a mathematically robust structure.

**4. Target Irreplaceability (10/10)**
I built a network model of all the RNA processing machinery in the parasite (17 components). Using spectral graph theory (the Fiedler vector), I calculated what happens when you remove each component. Removing the SL RNA causes a **62.6% drop in network connectivity** — the largest of any node. It's the single most critical component.

**5. Bayesian Design Optimization (6/10)**
I ran 201 Bayesian optimization evaluations to search for better ASO designs. Honest result: MRL-ASO-001 is **not** Pareto-optimal (rank 197/251). The optimal designs are all 30 nt. But at 25 nt, it's the best — and the length trade-off is justified by delivery considerations.

**6. Resistance Barrier (10/10)**
Using Markov chain models, I simulated how mutations could evolve in the SL RNA to escape the ASO. Result: **zero escape mutations** out of 75 analyzed positions. The worst-case time to resistance: **285 years**. The positions where the ASO binds are so critical that any mutation there kills the parasite before it can spread.

---

## Delivery: How Do You Get a Molecule Inside a Parasite Inside a Cell?

Here's the biological challenge: *Leishmania* parasites live **inside macrophages** (immune cells) of the infected dog. Getting MRL-ASO-001 to its target means crossing three barriers: the macrophage membrane, the phagolysosome membrane (an acidic compartment at pH 4.5), and the parasite membrane.

I built 6 computational modules to model each step of this journey:

**The Delivery Certificate: 60/60 (perfect score)**

- **Stability at pH 4.5**: The molecule remains 100% functional in the acidic phagolysosome, with a half-life of 1,083 hours (45 days). The LNA modifications keep the sugar puckers locked in the right conformation.

- **Membrane permeation**: Predicted intracellular concentration of 200,246 nM — over 2,000 times the minimum therapeutic threshold. Macrophages have an 18.7x uptake advantage over other cells (which is exactly what we want).

- **Conjugate strategy**: Attaching a trimannose sugar cluster targets the MRC1 receptor on macrophages, boosting uptake 9.7x. Macrophages *eat sugar* — we're exploiting their natural appetite.

- **LNP formulation**: An alternative delivery route using lipid nanoparticles (87 nm diameter, 98.8% encapsulation, 100% drug release at pH 4.5).

- **ADMET profile**: Therapeutic index of 8.0 (safe), 87% bioavailability, 21-day half-life (weekly dosing is feasible), preferential distribution to liver and spleen (where the parasites concentrate).

- **Immune simulation**: 1,000 Monte Carlo simulations of infected macrophages. Clearance probability: 100%. Mean time to clear parasites: 14 hours. EC50: 0.1 micromolar.

---

## Three Vaccine Platforms: From $81 to $0.00008 Per Dose

In parallel with the ASO work, I designed three vaccine platforms using the same 11 epitopes:

| Platform | System | Cost per dose | vs Leish-Tec |
|----------|--------|--------------|-------------|
| **mRNA-LNP** | In vivo expression (like COVID vaccines) | $5.75 | -85.8% |
| **E. coli** | Recombinant protein (pET-28a) | $2.28 | -91.6% |
| **L. tarentolae** | Live attenuated BSL-1 vector | $0.00008 | -100% |

The E. coli platform is the fastest path to market (24-36 months) and competes directly with Leish-Tec at 1/36th the cost.

The L. tarentolae platform is the cheapest vaccine I've ever seen modeled — less than a hundredth of a cent per dose. But there's a catch I documented honestly: the ASO designed to kill L. infantum would also kill the L. tarentolae vaccine vector (they share the same SL RNA in the binding region). Combined therapy requires a 2-week separation.

With 35.2 million dogs at risk in Brazil alone, the cost difference is the difference between mass vaccination and leaving most dogs unprotected.

---

## Track 3: When AI Enters the Lab

Once the foundation was solid, I built 11 AI/ML modules that analyze, extend, and connect everything:

### ESM-2: The Language of Proteins

I ran all 25 protein sequences through **ESM-2**, Meta's protein language model trained on 250 million proteins. It converts each protein into a vector of 1,280 numbers that capture deep evolutionary and structural properties.

Result: the model automatically clustered all 11 epitopes together and separated them from controls — confirming that our epitopes share something fundamentally in common that goes beyond what manual analysis can see.

### Reinforcement Learning: Teaching an AI to Design Better Molecules

I built an RL agent (REINFORCE algorithm) that learns to modify the ASO through trial and error:

- **State**: the current ASO sequence + chemical modifications
- **Actions**: change a base, toggle an LNA modification, switch backbone chemistry
- **Reward**: binding energy + melting temperature + GC content balance

After 500 training episodes (4 seconds on my laptop), the agent found variants with **28.5% higher reward** than the original MRL-ASO-001. The best variant has a binding energy of -36.83 kcal/mol.

But the AI Scientist module (more on that below) flagged a concern: all top RL variants have 60% GC content, which increases off-target risk. This is the kind of nuance that a pure optimization approach misses — and exactly why having a system that cross-checks its own results matters.

### Discrete Diffusion: Generating Novel Molecules

Inspired by Stable Diffusion (the image generator), I built a discrete diffusion model for biological sequences. Instead of generating images from noise, it generates molecules from noise.

The model produced 32 valid ASO variants with more conservative changes — the best has a binding energy of -29.25 kcal/mol with only 3 nucleotide changes from the original. These are closer to the original design but still novel.

### Sparse Autoencoders: Opening the Black Box

One of the most surprising results came from the Sparse Autoencoder (SAE). This module compresses protein information and forces the neural network to use very few active neurons — each one corresponding to an interpretable biological feature.

The top finding: **7 out of 10 most important features involve hydrophobicity**. N-terminal hydrophobicity emerged as the strongest predictor of whether a peptide is antigenic in Leishmania, with selectivity ratios up to 46x between epitopes and random controls.

This is a publishable finding that suggests a design principle: when engineering epitopes for Leishmania vaccines, prioritize hydrophobic residues at the N-terminus.

### The Digital Twin: Simulating a Dog

The Digital Twin module integrates three simulation layers into one:

1. **Pharmacokinetics**: how the drug distributes across organs (plasma, liver, spleen, kidney), accounting for disease-induced splenomegaly
2. **Immune dynamics**: how macrophages, parasites, IFN-gamma, and nitric oxide interact over time
3. **Stochastic variability**: each dog is different, so 100 virtual dogs are simulated with random noise

The result: a 28-day, hour-by-hour simulation of what happens when you inject MRL-ASO-001 into an infected dog.

### The AI Scientist: 5 Agents, 5 Hypotheses

The crown module reads all results from all other modules and runs 5 specialized agents:

- **LiteratureAgent**: searches 288 PubMed papers for relevant context
- **DesignAgent**: analyzes RL + diffusion outputs, proposes next designs
- **ValidationAgent**: cross-checks thermodynamic and structural consistency
- **KnowledgeAgent**: explores the Knowledge Graph for unexpected connections
- **ReportAgent**: synthesizes everything into a scientific narrative (powered by Claude)

Its latest output: **5 prioritized hypotheses** and **5 proposed experiments**, with the top recommendation being an RT-qPCR assay in L. infantum promastigotes — the simplest and most direct test of whether the ASO actually blocks trans-splicing.

---

## What I Didn't Find (Honest Negatives)

A platform that only reports successes is useless. Here are the honest failures:

- **MRL-003**: My custom small molecule bound human enzymes better than parasite enzymes. 20 redesigns failed. This is what killed the small molecule approach and led to the RNA pivot.

- **Oral drugs vs SL RNA**: I tested 15 oral compounds against the 3D structure of SL RNA. None bound. ASO is the only viable approach for this target.

- **Pareto optimality**: MRL-ASO-001 is NOT the mathematically optimal design. It ranks 197/251 in Bayesian optimization. The "better" designs are all 30 nt — longer, harder to deliver, more expensive to synthesize.

- **Sub-additive synergy**: I hypothesized that the dual mechanism (antisense + TLR9) would be synergistic. The simulation shows it's sub-additive — both mechanisms work, but 1+1 = 1.7, not 2+.

- **L. tarentolae orthogonality**: The cheapest vaccine platform has a fatal interaction with the ASO therapy. You can't give both simultaneously.

These aren't failures — they're data points that redirect the work.

---

## The Numbers

| Metric | Value |
|--------|-------|
| Total Python code | 84,000+ lines |
| Automated tests | 455 (all passing) |
| Result files | 80+ JSON files |
| Math validation | 52/60 VALIDATED |
| Delivery validation | 60/60 VALIDATED |
| AI/ML modules | 11 functional |
| PubMed papers indexed | 288 |
| ASO binding energy | -27.97 kcal/mol |
| Off-targets (human + dog) | 0 |
| Resistance timeline | 285 years |
| Host distance ratio | 14.5x |
| Vaccine cost (E. coli) | $2.28/dose |
| Market (Brazil) | R$880.75M/year |
| Time to build | Solo developer |

---

## What's Next

Everything described above is computational. Every single number is a prediction, not a measurement. The molecule hasn't been synthesized. No cell has been treated. No dog has been dosed.

The next step is the lab:

1. **Synthesize MRL-ASO-001** (~$500 from a commercial oligo synthesis company)
2. **RT-qPCR in promastigotes** (does the ASO actually block trans-splicing?)
3. **Binding assay (ITC/SPR)** (is the predicted -27.97 kcal/mol correct?)
4. **Canine macrophage assay** (does it kill parasites inside cells?)
5. **Mouse model** (does it work in a living organism?)

I'm currently enrolling in the PPGIT master's program at UFMG (Federal University of Minas Gerais) to pursue experimental validation, in collaboration with Prof. Dr. Rodolfo Giunchetti (ICB/UFMG), Dra. Wanessa Goes (CTVacinas), and Dra. Silvane Murta (Fiocruz).

---

## Why This Matters

Every year, thousands of dogs in Brazil are diagnosed with visceral leishmaniasis and euthanized because treatment options are inadequate. The disease also affects humans — and dogs are the primary reservoir.

A $2.28 vaccine and a targeted therapeutic with zero predicted off-targets won't solve this alone. Biology is messier than code. Predictions fail. Molecules that look perfect on screen can be useless in cells.

But the computational foundation is complete. The hypotheses are testable. The path from keyboard to lab bench is mapped.

Marley didn't make it. But maybe the next dog will.

---

*The Marley project is open-source and available at github.com/pedrohnsc2/marley. Built with Python, PyTorch, ESM-2, Claude, and an unreasonable amount of grief channeled into code.*

*If you're a researcher working on leishmaniasis, neglected tropical diseases, or antisense therapeutics and want to collaborate, reach out.*

---

**Pedro Nascimento** is a software developer and independent bioinformatics researcher based in Minas Gerais, Brazil. He builds computational tools for drug discovery and is pursuing a master's degree at PPGIT/UFMG focused on the experimental validation of the Marley platform.
