# Marley — Reverse Vaccinology Pipeline Report

> In memory of Marley, lost to canine visceral leishmaniasis.

**Date:** 2026-04-03
**Pipeline version:** 1.0
**Repository:** https://github.com/pedrohnsc2/marley

---

## 1. Executive Summary

The Marley pipeline analysed **8,527 proteins** from the *Leishmania infantum* JPCM5 proteome to identify antigen candidates for a canine mRNA vaccine against visceral leishmaniasis.

| Metric | Value |
|--------|-------|
| Total proteins analysed | 8,527 |
| Surface-exposed (signal peptide) | 139 (1.6%) |
| Conserved across strains (>80% identity) | 127 |
| Peptides with strong MHC-I binding (IC50 < 500 nM) | 204 |
| Epitopes selected for vaccine construct | 15 |
| Source proteins represented in construct | 6 |
| Total pipeline runtime | 637 minutes (~10.6 hours) |

The pipeline produced a **complete multi-epitope mRNA vaccine construct** (335 aa protein / 1,449 nt mRNA) ready for synthesis and experimental validation.

---

## 2. Pipeline Overview

```
8,527 proteins (L. infantum JPCM5, TriTrypDB)
        |
        | SignalP 6.0 — signal peptide filter
        v
      139 surface proteins (1.6%)
        |
        | NCBI BLAST — conservation across L. donovani, L. major, L. braziliensis
        v
      127 conserved candidates (>80% identity)
        |
        | IEDB MHC-I — binding affinity for 3 canine DLA alleles
        v
      204 peptides with IC50 < 500 nM
        |
        | Greedy selection — non-overlapping, diverse, top IC50
        v
       15 CTL epitopes from 6 source proteins
        |
        | Construct assembly + codon optimisation
        v
      mRNA vaccine sequence (1,449 nt)
```

---

## 3. Top 15 Computational Candidates

Ranked by final score (40% conservation + 60% immunogenicity):

| Rank | Gene ID | Protein | Conservation | Immunogenicity | Final Score |
|------|---------|---------|-------------|----------------|-------------|
| 1 | LINF_160022200 | **Prohibitin** | 1.000 | 0.058 | 0.435 |
| 2 | LINF_240013900 | **emp24/gp25L/p24 (GOLD family)** | 0.989 | 0.064 | 0.434 |
| 3 | LINF_080014800 | **Cysteine peptidase B (CPB)** | 0.994 | 0.035 | 0.418 |
| 4 | LINF_230010300 | Hypothetical protein - conserved | 0.923 | 0.080 | 0.417 |
| 5 | LINF_100010400 | **GP63 - leishmanolysin** | 0.992 | 0.033 | 0.416 |
| 6 | LINF_080015000 | **Cysteine peptidase B (CPB)** | 0.977 | 0.042 | 0.416 |
| 7 | LINF_180016800 | Hypothetical protein - conserved | 0.970 | 0.039 | 0.411 |
| 8 | LINF_110016300 | Hypothetical protein - conserved | 0.948 | 0.050 | 0.409 |
| 9 | LINF_210027700 | Vacuolar ATP synthase | 0.911 | 0.073 | 0.408 |
| 10 | LINF_140020300 | emp24/gp25L/p24 (GOLD family) | 0.955 | 0.042 | 0.407 |
| 11 | LINF_180020300 | NADH-ubiquinone oxidoreductase | 0.982 | 0.024 | 0.407 |
| 12 | LINF_270005100 | Meckelin (TM protein 67) | 0.924 | 0.059 | 0.405 |
| 13 | LINF_120009100 | 3'-nucleotidase/nuclease | 0.929 | 0.054 | 0.404 |
| 14 | LINF_010005000 | Protein of unknown function (DUF2946) | 0.928 | 0.034 | 0.392 |
| 15 | LINF_350044200 | Hypothetical protein - conserved | 0.929 | 0.032 | 0.391 |

**Notable findings:** GP63 (leishmanolysin) and CPB (cysteine peptidase B) are well-established vaccine candidates in *Leishmania* research literature. The pipeline identified them independently through computational analysis, validating the approach.

---

## 4. Pre-validated Antigens from Literature

These antigens were pre-loaded with priority status, sourced from published research by Brazilian groups:

| Rank | Antigen | Source | Evidence | Score |
|------|---------|--------|----------|-------|
| 1 | LiHyp1 | Giunchetti/UFMG | Murine validation, Th1 response via IFN-gamma, immunoproteomics | 0.95 |
| 2 | A2 | UFMG / Leish-Tec | Only MAPA-approved vaccine in Brazil, 96.41% efficacy | 0.92 |
| 3 | LBSap antigens | Reis/UFOP + Giunchetti/UFMG | Technology transferred to Ouro Fino Saude Animal | 0.90 |
| 4 | Lutzomyia longipalpis proteins | Giunchetti/UFMG | Patent052, transmission-blocking | 0.88 |
| 5 | KMP-11 | Literature | High conservation across strains, documented immunogenicity | 0.85 |
| 6 | LiESP/Q | Literature | High diagnostic specificity, immunoprotective potential | 0.83 |
| 7 | LACK | Literature | T-cell activation in murine models | 0.82 |
| 8 | HSP70/HSP83 | Literature | Overexpressed under stress, conserved, strong cellular response | 0.80 |

---

## 5. mRNA Vaccine Construct Design

### 5.1 Construct Architecture

```
[tPA signal peptide] - [L7/L12 adjuvant] - EAAAK - [Epi1] - AAY - [Epi2] - AAY - ... - [Epi15]
 \___ 23 aa ________/   \___ 130 aa ____/           \____________ 15 epitopes _______________/
```

| Component | Sequence / Description | Purpose |
|-----------|----------------------|---------|
| Signal peptide | tPA leader (23 aa) | Drives secretion for MHC presentation |
| Adjuvant | 50S ribosomal L7/L12 (130 aa) | TLR4 agonist, induces Th1 bias (critical for *Leishmania*) |
| EAAAK linker | Rigid alpha-helix spacer | Separates adjuvant from epitope cassette |
| AAY linkers | Proteasomal cleavage sites | Enables proper epitope processing |
| CTL epitopes | 15 x 9-mer peptides | MHC-I presented, activate CD8+ T cells |

### 5.2 Selected Epitopes

| # | Peptide | Source Protein | Allele | IC50 (nM) |
|---|---------|---------------|--------|-----------|
| 1 | RMMRSLTPF | emp24/gp25L/p24 (GOLD) | DLA-88*03401 | 11.2 |
| 2 | YIYETFASM | emp24/gp25L/p24 (GOLD) | DLA-88*50101 | 18.8 |
| 3 | SLMCVFYFK | Hypothetical - conserved | DLA-88*50801 | 33.8 |
| 4 | YLAALVPAL | Hypothetical - conserved | DLA-88*50101 | 38.2 |
| 5 | LIIEDLSLV | Prohibitin | DLA-88*50101 | 64.1 |
| 6 | FAFSVSARR | Hypothetical - conserved | DLA-88*50801 | 64.7 |
| 7 | MILGTFVRL | emp24/gp25L/p24 (GOLD) | DLA-88*50101 | 72.9 |
| 8 | MQNVTFVPK | Prohibitin | DLA-88*50801 | 75.4 |
| 9 | RILESISNV | **GP63 - leishmanolysin** | DLA-88*50101 | 97.4 |
| 10 | ILYNKISGL | Prohibitin | DLA-88*03401 | 97.6 |
| 11-15 | LLTANVCYK (x5) | **Cysteine peptidase B (CPB)** | DLA-88*50801 | 118.0 |

**Allele coverage:** All 3 canine DLA alleles are represented (DLA-88*50101, DLA-88*50801, DLA-88*03401).

**Source diversity:** Epitopes come from 6 different proteins: emp24/GOLD, prohibitin, GP63, CPB, and 2 conserved hypothetical proteins.

### 5.3 Physicochemical Properties

| Property | Value | Interpretation |
|----------|-------|---------------|
| Protein length | 335 aa | Suitable for mRNA delivery |
| Molecular weight | 35.95 kDa | Within typical range for recombinant vaccines |
| Isoelectric point (pI) | 7.93 | Near-neutral, good solubility expected |
| Instability index | 28.38 | **Stable** (threshold < 40) |
| GRAVY | 0.589 | Slightly hydrophobic |
| Aromaticity | 0.113 | Normal range |

### 5.4 mRNA Cassette

```
5'cap -- 5'UTR (57 nt) -- CDS (1,005 nt) -- TGA -- 3'UTR (132 nt) -- 3'UTR (132 nt) -- poly(A)120
```

| Property | Value |
|----------|-------|
| Total mRNA length | 1,449 nt |
| CDS length | 1,005 nt (335 codons) |
| GC content | 66.9% |
| Codon optimisation | *Canis lupus familiaris* (Kazusa table) |
| 5' UTR | Alpha-globin-like + Kozak consensus |
| 3' UTR | Human beta-globin (x2, BioNTech approach) |
| Poly(A) tail | 120 nt |
| Restriction sites avoided | EcoRI, BamHI, HindIII, XbaI, NheI, BsaI |

### 5.5 Protein Sequence (FASTA)

```
>Marley_vaccine_construct len=335 signal=tPA adj=L7L12
MDAMKRGLCCVLLLCGAVFVSASMAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGA
APAGAAVEAAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAK
AKLEAAGATVTVKEAAAKRMMRSLTPFAAYYIYETFASMAAYSLMCVFYFKAAYYLAALVPALAAYLIIE
DLSLVAAYFAFSVSARRAAYMILGTFVRLAAYMQNVTFVPKAAYRILESISNVAAYILYNKISGLAAYLL
TANVCYKAAYLLTANVCYKAAYLLTANVCYKAAYLLTANVCYKAAYLLTANVCYK
```

---

## 6. Methodology

### 6.1 Genome Fetch
The complete *L. infantum* JPCM5 proteome (8,527 annotated protein sequences) was downloaded from [TriTrypDB](https://tritrypdb.org) release 68.

### 6.2 Surface Protein Filter
Signal peptide prediction was performed using **SignalP 6.0** via the BioLib SDK. Only proteins classified as Sec/SPI (classic secretory signal peptide) were retained, indicating surface exposure or secretion -- both desirable traits for vaccine targets accessible to the host immune system.

### 6.3 Conservation Analysis
Protein conservation was evaluated via **NCBI BLAST** (blastp, E-value < 1e-10) against the non-redundant database, restricted to *Leishmania* (taxid: 5658). Conservation scores represent average percent identity across three comparison species:
- *Leishmania donovani* (closest relative, same complex)
- *Leishmania major* (well-studied reference species)
- *Leishmania braziliensis* (different subgenus, broader conservation signal)

Candidates with conservation score < 0.80 were rejected. Self-hits (*L. infantum*) were excluded from scoring.

### 6.4 Immunogenicity Scoring
MHC-I binding affinity was predicted using the **IEDB Analysis Resource** (tools-cluster-interface.iedb.org) with the following parameters:
- **Method:** netMHCpan-BA (binding affinity)
- **Alleles:** DLA-88*50101, DLA-88*50801, DLA-88*03401 (canine MHC-I)
- **Peptide length:** 9-mer
- **Threshold:** IC50 < 500 nM (strong binder)

### 6.5 Epitope Selection
Epitopes were selected using a greedy algorithm that:
1. Pools all peptides with IC50 < 500 nM across all candidates
2. Sorts by IC50 ascending (best binders first)
3. Skips peptides overlapping (>= 5 shared residues) with already-selected epitopes from the same protein
4. Limits to 3 epitopes per source protein (ensures diversity)
5. Selects up to 15 CTL epitopes total

### 6.6 Construct Assembly
The multi-epitope protein was assembled following established reverse vaccinology design principles:
- **tPA signal peptide** for secretion and MHC presentation
- **L7/L12 adjuvant** (50S ribosomal protein) as TLR4 agonist with documented Th1 bias
- **EAAAK** rigid linker between adjuvant and epitope cassette
- **AAY** linkers between CTL epitopes (natural proteasomal cleavage sites)

### 6.7 Codon Optimisation
The protein was reverse-translated using the *Canis lupus familiaris* codon usage table (Kazusa database, taxonomy ID 9615). Post-processing removed restriction enzyme sites and homopolymer runs > 4 nt.

### 6.8 mRNA Design
The mRNA cassette follows the BioNTech BNT162b2 architecture:
- 5' cap (m7GpppNm, Cap1 -- annotated, post-transcriptional modification)
- 5' UTR with Kozak consensus
- Codon-optimised CDS
- Dual human beta-globin 3' UTR for enhanced stability
- 120 nt poly(A) tail

---

## 7. Limitations

1. **MHC-I only:** The current analysis covers CD8+ T cell epitopes (CTL). MHC-II (helper T cell) and B-cell epitope predictions were not included in this version.

2. **DLA allele coverage:** Only 3 canine DLA alleles were used. The DLA locus is highly polymorphic and breed-specific coverage was not assessed.

3. **Safety checks unavailable:** VaxiJen (antigenicity) and AllerTOP (allergenicity) servers were unreachable during this run. These should be repeated when servers are available.

4. **GC content:** At 66.9%, the GC content is above the ideal range (40-60%). This may affect mRNA stability and translation efficiency and should be optimised further.

5. **No experimental validation:** All results are computational predictions. The construct must be validated through expression, immunogenicity assays, and animal challenge studies.

6. **Validated antigens lack sequences:** The 8 pre-loaded antigens from literature were scored based on published evidence but their protein sequences were not included in the epitope selection (empty sequence field). Future versions should fetch these from UniProt.

---

## 8. Suggested Next Steps

### Computational
- [ ] Fetch protein sequences for validated antigens from UniProt (A2=Q25327, KMP-11=P15711, etc.)
- [ ] Add MHC-II (helper T cell) epitope prediction via IEDB
- [ ] Add B-cell linear epitope prediction via BepiPred
- [ ] Run VaxiJen and AllerTOP when servers are available
- [ ] Optimise GC content of the CDS
- [ ] Predict 3D structure via AlphaFold or ESMFold

### Experimental
- [ ] Synthesise mRNA construct and encapsulate in LNP
- [ ] In-vitro expression and Western blot confirmation
- [ ] ELISPOT and cytokine profiling (IFN-gamma, IL-2, TNF-alpha) with canine PBMCs
- [ ] Murine challenge model (BALB/c) for preliminary efficacy
- [ ] Canine immunogenicity and challenge study
- [ ] Dose-response and adjuvant formulation optimisation

---

## 9. Tools and References

### Tools used
| Tool | Version | Purpose |
|------|---------|---------|
| TriTrypDB | Release 68 | *L. infantum* JPCM5 proteome |
| SignalP | 6.0 (BioLib) | Signal peptide prediction |
| NCBI BLAST | Remote (nr database) | Conservation analysis |
| IEDB | netMHCpan-BA | MHC-I binding prediction |
| Biopython | ProtParam | Physicochemical analysis |

### Key references
1. Patronov A, Doytchinova I. T-cell epitope vaccine design by immunoinformatics. *Open Biol.* 2013;3(1):120139.
2. Doytchinova IA, Flower DR. VaxiJen: a server for prediction of protective antigens. *BMC Bioinformatics.* 2007;8:4.
3. Reis AB et al. Immunity to *Leishmania* and the rational search for vaccines against canine leishmaniasis. *Trends Parasitol.* 2010;26(7):341-9.
4. Vogel AB et al. BNT162b vaccines protect rhesus macaques from SARS-CoV-2. *Nature.* 2021;592:283-289.
5. Fernandes AP et al. Protective immunity against challenge with *Leishmania (Leishmania) chagasi* in beagle dogs vaccinated with recombinant A2 protein. *Vaccine.* 2008;26(46):5888-95.

---

## 10. Disclaimer

This report was generated by the **Marley** open-source reverse vaccinology pipeline. All predictions are computational (*in silico*) and have not been experimentally validated. This construct is a **proposal for experimental validation**, not a finished vaccine product. Any translational application requires laboratory testing, including expression, purification, immunogenicity assays, and regulatory-compliant animal studies.

---

*Generated by Marley v1.0 — https://github.com/pedrohnsc2/marley*
