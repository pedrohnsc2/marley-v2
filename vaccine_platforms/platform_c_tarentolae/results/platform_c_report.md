# Platform C -- L. tarentolae Live Vaccine (LEXSY)

**Generated:** 2026-04-11T00:26:02.824625+00:00
**Pipeline version:** 0.1.0
**Orthogonality check:** [WARNING] WARNING

---

## 0. MANDATORY: SL RNA Orthogonality Check [WARNING]

**Purpose:** Verify that ASO MRL-ASO-001 (designed to kill L. infantum) does
not also kill L. tarentolae (our expression vector).

### Sequence Comparison

| Property | L. infantum | L. tarentolae |
|---|---|---|
| SL RNA (39 nt) | `AACTAACGCTATATAAGTATCAGTTTCTGTACTTATATG` | `AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG` |
| ASO target region | `ACGCTATATAAGTATCAGTTTCTGT` | `ACGCTATATAAGTATCAGTTTCTGT` |

### Alignment

| Metric | Value |
|---|---|
| Alignment length | 25 nt |
| Matches (ASO vs L. tarentolae) | 25 |
| Mismatches | 0 |
| Mismatch positions | [] |
| Complementarity to L. infantum | 100.0% |
| Complementarity to L. tarentolae | 100.0% |

### Thermodynamics (SantaLucia 1998)

| Parameter | L. infantum | L. tarentolae | Difference |
|---|---|---|---|
| DeltaH (kcal/mol) | -119.90 | -119.90 | -- |
| DeltaS (cal/mol/K) | -338.80 | -338.80 | -- |
| Tm (C) | 49.37 | 49.37 | 0.00 |
| DeltaG at 37C (kcal/mol) | -14.82 | -14.82 | 0.00 |

**Selectivity:** DeltaDeltaG = 0.00 kcal/mol (threshold: >= 2.0 kcal/mol)
**Kd ratio:** 1.0x (L. tarentolae Kd / L. infantum Kd)

### Verdict

> WARNING -- Diferenca de DeltaG insuficiente (0.00 kcal/mol < 2.0 kcal/mol). O ASO MRL-ASO-001 pode afetar a viabilidade de L. tarentolae como vetor de expressao. REQUER AVALIACAO EXPERIMENTAL.

**Warnings:**
- CRITICO: DeltaDeltaG = 0.00 kcal/mol esta abaixo do limiar de 2.0 kcal/mol. Isso significa que o ASO pode ter afinidade significativa pelo SL RNA de L. tarentolae, potencialmente inibindo o trans-splicing do vetor de expressao.
- ALERTA: Zero mismatches detectados -- o ASO tem complementaridade perfeita com ambos os SL RNAs. A seletividade depende inteiramente de fatores cineticos e de acessibilidade estrutural.

---

## 1. Why L. tarentolae? (C1)

| Aspect | Detail |
|---|---|
| Biosafety | BSL-1 (nao patogenico para mamiferos; infecta lagartixas) |
| Growth conditions | 26C, meio BHI suplementado, sem CO2, tempo de duplicacao ~8h |
| PTM capability | N-glicosilacao, ancoras GPI, acetilacao N-terminal -- equivalentes a L. infantum patogenica |
| Intrinsic adjuvant | Moleculas do parasita (LPG, GP63, GIPLs) estimulam TLR2/TLR4 do sistema imune inato, eliminando necessidade de adjuvante exogeno para a forma viva |
| Literature precedent | HPV (Breitling 2002), HCV (Niimi 2018), SARS-CoV-2 (Pirdel 2021) -- proteinas vacinais expressas com sucesso em L. tarentolae |
| Commercial system | LEXSY (Jena Bioscience): kit comercial maduro com vetores, cepas e protocolos otimizados para expressao constitutiva e induzivel |
| Vector | pLEXSY-sat2 |
| Selection | nourseothricin (100 ug/mL) |

---

## 2. Construct Redesign (C2)

### Architecture

```
[SAP1 signal] - [L7/L12 adjuvant] - EAAAK - [11 epitopes + AAY] - GSGSGS - [His6] - [Strep-tag II]
```

Genomic integration context:
```
[SSU 5' flank] - [5'UTR] - ATG[CDS]TGA - [3'UTR] - [SSU 3' flank]
```

| Property | Value |
|---|---|
| Total protein length | 302 aa |
| Molecular weight | 32,580.50 Da (32.6 kDa) |
| Isoelectric point (pI) | 6.21 |
| Instability index | 30.58 (stable) |
| GRAVY | 0.4086 |
| Aromaticity | 0.1159 |
| Unique epitopes | 11 |
| Vector | pLEXSY-sat2 |
| Selection | nourseothricin |

### Modifications from mRNA Platform

| Component | mRNA (Platform A) | L. tarentolae (Platform C) |
|---|---|---|
| signal_peptide | tPA (MDAMKRGLCCVLLLCGAVFVSAS) | SAP1 (MKFFVFALFAVALCSAEA) |
| purification_tag | None (mRNA traduzido in vivo) | His6-GSGSGS-Strep (HHHHHHGSGSGSWSHPQFEK) |
| utrs | 5'UTR humano + 3'UTR com poli(A) | Regioes intergenicas de L. tarentolae para trans-s |
| integration | None (mRNA transiente) | SSU rRNA locus (recombinacao homologa) |
| codon_usage | Canis lupus familiaris (otimizado para traducao in | L. tarentolae (GC-rich, trypanosomatid bias) |
| vector | None (mRNA sintetico) | pLEXSY-sat2 (nourseothricin) |

### Epitope Cassette

| # | Peptide | Linker |
|---|---|---|
| 1 | `RMMRSLTPF` | AAY |
| 2 | `YIYETFASM` | AAY |
| 3 | `SLMCVFYFK` | AAY |
| 4 | `YLAALVPAL` | AAY |
| 5 | `LIIEDLSLV` | AAY |
| 6 | `FAFSVSARR` | AAY |
| 7 | `MILGTFVRL` | AAY |
| 8 | `MQNVTFVPK` | AAY |
| 9 | `RILESISNV` | AAY |
| 10 | `ILYNKISGL` | AAY |
| 11 | `LLTANVCYK` | (C-terminal -> tags) |

### Protein Sequence

```
MKFFVFALFAVALCSAEAMAKLSTDELLDAFKEMTLLELSDFVKKFEETFEVTAAAPVAVAAAGAAPAGA
AVEAAEEQSEFDVILEAAGDKKIGVIKVVREIVSGLGLKEAKDLVDGAPKPLLEKVAKEAADEAKAKLEA
AGATVTVKEAAAKRMMRSLTPFAAYYIYETFASMAAYSLMCVFYFKAAYYLAALVPALAAYLIIEDLSLV
AAYFAFSVSARRAAYMILGTFVRLAAYMQNVTFVPKAAYRILESISNVAAYILYNKISGLAAYLLTANVC
YKGSGSGSHHHHHHWSHPQFEK
```

---

## 3. Codon Optimization for L. tarentolae (C3)

| Metric | Value | Target |
|---|---|---|
| DNA insert length | 906 nt | -- |
| GC content | 67.8% | 55%-65% |
| GC in range | NO | Yes |
| CAI | 0.9992 | >0.85 |
| CAI above threshold | Yes | Yes |
| Rare codons (< 10%) | 0 | 0 |
| Restriction sites | 0 | 0 |
| Homopolymer runs | 0 | 0 |
| Poly-A signal | None | None |
| Poly-T signal | None | None |

**Warnings:**
- GC content = 67.8% fora da faixa alvo (55%-65%)

### DNA Insert Sequence (CDS only)

```
ATGAAGTTCTTCGTGTTCGCCCTGTTCGCCGTGGCCCTGTGCAGCGCCGAGGCGATGGCCAAGCTGAGCA
CCGACGAGCTGCTGGACGCCTTCAAGGAGATGACCCTGCTGGAGCTGAGCGACTTCGTGAAGAAGTTCGA
GGAGACCTTCGAGGTGACCGCCGCCGCCCCGGTGGCCGTGGCCGCCGCCGGCGCCGCCCCGGCCGGCGCC
GCCGTGGAGGCCGCCGAGGAGCAGAGCGAGTTCGACGTGATCCTGGAGGCCGCCGGCGACAAGAAGATCG
GCGTGATCAAGGTGGTGCGCGAGATCGTGAGCGGCCTGGGCCTGAAGGAGGCCAAGGACCTGGTGGACGG
CGCCCCGAAGCCGCTGCTGGAGAAGGTGGCCAAGGAGGCCGCCGACGAGGCCAAGGCCAAGCTGGAGGCC
GCCGGCGCCACCGTGACCGTGAAGGAGGCCGCCGCCAAGCGCATGATGCGCAGCCTGACCCCGTTCGCCG
CCTACTACATCTACGAGACCTTCGCCAGCATGGCCGCCTACAGCCTGATGTGCGTGTTCTACTTCAAGGC
CGCCTACTACCTGGCCGCCCTGGTGCCGGCCCTGGCCGCCTACCTGATCATCGAGGACCTGAGCCTGGTG
GCCGCCTACTTCGCCTTCAGCGTGAGCGCCCGCCGCGCCGCCTACATGATCCTGGGCACCTTCGTGCGCC
TGGCCGCCTACATGCAGAACGTGACCTTCGTGCCGAAGGCCGCCTACCGCATCCTGGAGAGCATCAGCAA
CGTGGCCGCCTACATCCTGTACAACAAGATCAGCGGCCTGGCCGCCTACCTGCTGACCGCCAACGTGTGC
TACAAGGGCAGCGGCAGCGGCAGCCACCACCACCACCACCACTGGAGCCACCCGCAGTTCGAGAAG
```

---

## 4. Cost Model (C4)

### Modalidade A: Vacina Viva (Organismo Inteiro)

| Item | Value |
|---|---|
| Growth medium cost | $18.00/L |
| Doses per liter | 50,000 |
| **Cost per dose** | **$0.000660** |
| **Cost per animal (3 doses)** | **$0.001980** |
| Industrial cost per dose | $0.000079 |
| Adjuvant | intrinsic (parasite molecules) |
| Cold chain | 2-8C (cultura viva) ou liofilizado (temp ambiente) |

### Modalidade B: Proteina Secretada Purificada

| Item | Value |
|---|---|
| Growth + purification cost | $113.00/L |
| Yield (purified) | 6.0 mg/L |
| Cost per mg | $21.33 |
| Protein cost per dose | $1.0667 |
| **Total cost per dose** | **$3.0667** |
| **Cost per animal (3 doses)** | **$9.2000** |
| Industrial total per dose | $2.1280 |
| Adjuvant | QuilA ($1.50/dose) |
| Cold chain | 2-8C (proteina purificada) |

### Setup Costs (One-Time)

| Item | Cost |
|---|---|
| pLEXSY kit | $1,500.00 |
| Amortized per batch | $15.00 |
| BSL level | BSL-1 |
| Growth temperature | 26C |

### Platform Comparison

| Platform | Cost/dose | Cold chain | Adjuvant |
|---|---|---|---|
| platform_a_mrna | $5.00-15.00 | -20C a -80C | LNP (nanoparticula lipidica) incluida |
| platform_b_ecoli | $2.00-4.00 | 2-8C | QuilA ($1.50/dose) |
| platform_c_tarentolae_live | $0.0007 | 2-8C (cultura viva) ou liofilizado (temp ambiente) | Intrinseco (moleculas do parasita) |
| platform_c_tarentolae_protein | $3.0667 | 2-8C (proteina purificada) | QuilA ($1.50/dose) |

---

## 5. ASO Compatibility Assessment (C5)

### Can MRL-ASO-001 and L. tarentolae vaccine be used together?

**Answer:** COM RESTRICOES -- a diferenca de afinidade e pequena (DeltaDeltaG = 0.00 kcal/mol, abaixo do limiar de 2.0 kcal/mol). O ASO pode afetar a viabilidade de L. tarentolae.

### Recommended Protocol

**Protocolo com precaucoes:**

1. **NUNCA administrar simultaneamente.** O ASO pode matar L. tarentolae.

2. **Opcao A -- Vacinar primeiro:** Administrar vacina L. tarentolae, aguardar resposta imune (2-4 semanas), depois iniciar ASO.

3. **Opcao B -- Intervalo longo:** ASO primeiro, aguardar >1 semana de washout antes da vacinacao.

4. **Opcao C -- Usar Platform B (E. coli):** Se a compatibilidade com ASO e prioritaria, a vacina recombinante em E. coli nao depende de SL RNA e e totalmente compativel.

5. **Opcao D -- Modificar o ASO:** Introduzir modificacoes quimicas (LNA, 2'-OMe) nas posicoes de mismatch para aumentar seletividade por L. infantum.

### Risk Assessment

| Factor | Assessment |
|---|---|
| Thermodynamic selectivity | DeltaDeltaG = 0.00 kcal/mol (ABAIXO do limiar de 2.0) |
| Mismatch protection | Apenas 0 mismatch(es) -- protecao insuficiente |
| Temporal separation | Requer intervalo prolongado (>1 semana) entre ASO e vacina |
| Overall risk | MODERADO a ALTO -- requer precaucoes ou plataforma alternativa |

### Clinical Implications

> INCOMPATIVEL sem modificacoes. Opcoes: (1) Administrar vacina L. tarentolae ANTES do tratamento com ASO; (2) Usar intervalo longo (>1 semana) entre ASO e vacinacao; (3) Modificar o ASO com LNA/2'-OMe nas posicoes de mismatch para aumentar seletividade; (4) Considerar vetor de expressao alternativo (E. coli, Platform B).

---

*Report generated by the Marley reverse vaccinology pipeline -- Platform C (L. tarentolae LEXSY)*
