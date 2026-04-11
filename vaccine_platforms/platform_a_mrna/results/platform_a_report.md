# Platform A -- mRNA Vaccine Optimization

**Generated:** 2026-04-11T00:26:02.780759+00:00
**Pipeline version:** 0.1.0
**Construct:** 335 aa, 1449 nt

---

## 1. Epitope Redundancy Analysis (A1)

### Overview

| Metric | Value |
|---|---|
| Total canonical epitopes | 11 |
| Unique source genes | 5 |
| Unique DLA alleles | 3 |
| Pairwise comparisons | 55 |
| Redundant pairs (>70% sim + same allele) | 0 |
| Minimum epitopes for 95% coverage | 11 |

### Top 10 Most Similar Pairs

| Epitope A | Epitope B | Similarity | Same Gene | Same Allele |
|---|---|---|---|---|
| `LIIEDLSLV` | `RILESISNV` | 44.44% | No | Yes |
| `MILGTFVRL` | `MQNVTFVPK` | 44.44% | No | No |
| `YIYETFASM` | `MILGTFVRL` | 33.33% | Yes | Yes |
| `YLAALVPAL` | `LLTANVCYK` | 33.33% | No | No |
| `RMMRSLTPF` | `RILESISNV` | 22.22% | No | No |
| `YIYETFASM` | `LIIEDLSLV` | 22.22% | No | Yes |
| `YIYETFASM` | `MQNVTFVPK` | 22.22% | No | No |
| `YIYETFASM` | `RILESISNV` | 22.22% | No | Yes |
| `SLMCVFYFK` | `MQNVTFVPK` | 22.22% | No | Yes |
| `SLMCVFYFK` | `LLTANVCYK` | 22.22% | No | Yes |

### Construct Variants

| Epitopes | Length (aa) | MW (kDa) | II | Stable | Antigenicity | Coverage |
|---|---|---|---|---|---|---|
| 11 | 287 | 30.7 | 32.0 | Yes | 0.6245 | 100.0% |
| 9 | 263 | 28.1 | 32.7 | Yes | 0.6231 | 92.4% |
| 7 | 239 | 25.4 | 31.6 | Yes | 0.6270 | 84.1% |

**Recommendation:** Nenhum par redundante encontrado. Todos os 11 epitopos sao suficientemente distintos em sequencia e especificidade alelica. Recomenda-se manter o construto completo com 11 epitopos.

---

## 2. Pseudouridine Optimization (A2)

### mRNA Characteristics

| Metric | Value |
|---|---|
| mRNA length | 1008 nt |
| Total uridines | 154 |
| Uridine fraction | 15.3% |
| GC content | 67.3% |

### Uridine Structural Context

| Category | Count |
|---|---|
| High base-pair score (>0.5) | 5 |
| Medium base-pair score (0.25-0.5) | 149 |
| Low base-pair score (<0.25) | 0 |
| Average base-pair score | 0.4198 |

### Substitution Strategies

| Strategy | U Substituted | Translation Fold | Immunogenicity Reduction | Cost Factor |
|---|---|---|---|---|
| 30% psiU (selective top-ranked) | 46/154 | 1.65x | 25.5% | 1.15x |
| 50% psiU (selective top-ranked) | 77/154 | 1.92x | 42.5% | 1.25x |
| 100% psiU (full substitution) | 154/154 | 2.50x | 85.0% | 1.50x |

**Optimal fraction:** 50%

**Rationale:** A substituicao de 50% das uridinas (priorizando posicoes em stems de estrutura secundaria) oferece o melhor equilibrio entre: (1) aumento de ~1.75x na eficiencia de traducao, (2) reducao de ~42% na ativacao imune inata (TLR7/8), (3) custo incremental moderado de ~1.25x. A substituicao total (100%) maximiza traducao (~2.5x) e reduz imunogenicidade em ~85%, mas o custo de m1psi-UTP e significativo para producao em escala veterinaria. Para um primeiro candidato pre-clinico, 50% seletivo e a abordagem mais custo-efetiva. Na fase clinica, 100% pode ser justificado pelo beneficio imunologico.

---

## 3. Veterinary LNP Formulation (A3)

### Dose Calculation

| Parameter | Value |
|---|---|
| Target species | Canis lupus familiaris |
| Body weight | 15.0 kg |
| Dose per kg | 30.0 ug/kg |
| Total dose per injection | 450.0 ug |
| Doses per animal | 2 |
| Total mRNA per animal | 900.0 ug |

### LNP Composition (per dose)

| Component | Mol% | Amount |
|---|---|---|
| SM-102 (ionizable lipid) | 50.0% | 8181.73 nmol |
| DSPC | 10.0% | 1636.35 nmol |
| Cholesterol | 38.5% | 6299.94 nmol |
| PEG-DMG 2000 | 1.5% | 245.45 nmol |

| Metric | Value |
|---|---|
| N/P ratio | 6.0 |
| Total lipid mass | 9146.16 ug |
| Lipid:mRNA ratio (w/w) | 20.32 |
| Encapsulation efficiency | 95.0% |

### N/P Ratio Optimization

| N/P | Lipid Mass (ug) | Encapsulation | Transfection Score | Optimal |
|---|---|---|---|---|
| 6 | 9146.16 | 96.3% | 0.980 | **Yes** |
| 7 | 10670.52 | 96.7% | 0.995 |  |
| 8 | 12194.88 | 97.0% | 1.000 |  |
| 9 | 13719.23 | 96.7% | 0.995 |  |
| 10 | 15243.59 | 96.3% | 0.980 |  |
| 11 | 16767.95 | 96.0% | 0.855 |  |
| 12 | 18292.31 | 95.7% | 0.720 |  |

### Storage Options

| Condition | Temperature | Shelf Life | Feasibility |
|---|---|---|---|
| Ultra-frozen (-70 C) | -70.0 C | 6 months | Laboratorial e centros urbanos |
| Frozen (-20 C) | -20.0 C | 3 months | Ampla distribuicao urbana |
| Refrigerated (2-8 C, lyophilized) | 5.0 C | 12 months | Distribuicao ampla incluindo areas rurais |

### Cost Estimate

| Component | Laboratory | Industrial |
|---|---|---|
| IVT (mRNA production) | $225.00 | $2.2500 |
| Lipids | $0.0127 | $0.000635 |
| Formulation | $15.00 | $1.50 |
| QC | $20.00 | $2.00 |
| **Total per dose** | **$260.01** | **$5.7506** |
| **Total per animal** | **$520.03** | **$11.5013** |

**Recommendation:** Para um cao de 15 kg, recomenda-se 450 ug de mRNA por dose (30 ug/kg), formulado em LNP com composicao SM-102:DSPC:Chol:PEG-DMG (50:10:38.5:1.5 mol%) e razao N/P = 6. Regime de 2 doses (primo + reforco). Para distribuicao em areas endemicas brasileiras, a formulacao liofilizada (2-8 C, 12 meses de validade) e a mais adequada, apesar do custo adicional de liofilizacao. Custo estimado por animal: $520.03 (lab) / $11.50 (industrial).

---

## 4. Strategic Funding & Market Analysis (A4)

### Market Size (Brazil)

| Metric | Value |
|---|---|
| Total dogs in Brazil | 54,200,000 |
| Dogs at risk (endemic areas) | 35,230,000 (65.0%) |
| Mean seroprevalence | 25.0% |
| Estimated infected dogs | 8,807,500 |
| Addressable market | 7,046,000 dogs/year |
| Market value (BRL/year) | R$ 880,750,000 |
| Market value (USD/year) | $ 158,535,000 |
| Human cases/year | 3,500 (CFR: 8.0%) |
| Government control cost | R$ 350,000,000/year |

### Vaccine Sovereignty Argument

**Soberania vacinal em mRNA: do COVID-19 a One Health**

A pandemia de COVID-19 expôs a dependencia do Brasil de tecnologia estrangeira para vacinas de mRNA. O pais foi um dos ultimos a acessar vacinas de mRNA, com impacto direto em mortalidade. Desde entao, o governo brasileiro investiu em capacidade nacional de mRNA.

**COVID-19 precedent:** Em 2023, Bio-Manguinhos/Fiocruz firmou acordo de transferencia de tecnologia com BioNTech para producao de vacinas de mRNA no Brasil. A planta de producao em Santa Cruz (RJ) esta em fase de instalacao. O investimento total excede R$ 1 bilhao.

**One Health:** A leishmaniose visceral e uma zoonose: o cao e o principal reservatorio urbano de Leishmania infantum. Vacinar caes reduz a prevalencia canina e, consequentemente, a transmissao para humanos. O PNCLV (Programa Nacional de Controle da LV) do Ministerio da Saude atualmente depende de eutanasia de caes soropositivos -- estrategia eticamente controversa e de eficacia limitada. Uma vacina eficaz transformaria o programa.

### Funding Sources (ranked by score)

| Rank | Source | Score | Probability | Amount (BRL) | Timeline |
|---|---|---|---|---|---|
| 1 | Fiocruz / Biomanguinhos -- Parceria Tecnologica | 8.0/10 | 20% | R$ 5.0-50.0M | 18-36 months |
| 2 | FINEP Inovacao | 7.5/10 | 40% | R$ 0.5-10.0M | 8-18 months |
| 3 | MCTI / CNPq -- Chamada Universal + INCT | 7.0/10 | 55% | R$ 0.1-2.0M | 6-12 months |
| 4 | BNDES Saude Animal | 6.5/10 | 35% | R$ 2.0-20.0M | 12-24 months |
| 5 | FAPESP / FAPs estaduais -- PIPE/PAPPE | 6.0/10 | 45% | R$ 0.2-2.5M | 6-18 months |
| 6 | USDA SBIR / International Collaborations | 4.5/10 | 15% | R$ 0.8-5.0M | 12-24 months |

### Regulatory Pathways

**Brazil (MAPA):** 3-5 anos
**USA (USDA-CVB):** 2-5 anos

### Strategic Recommendation

Estrategia de financiamento em tres fases:

FASE 1 (0-12 meses) -- Prova de conceito: ~R$ 500K-2M
  Fonte primaria: MCTI / CNPq -- Chamada Universal + INCT (CNPq/MCTI)
  Fonte secundaria: BNDES Saude Animal (FAP estadual)
  Objetivo: dados pre-clinicos em camundongos BALB/c

FASE 2 (12-30 meses) -- Pre-clinico avancado: ~R$ 5-10M
  Fonte primaria: FINEP Inovacao (FINEP subvencao)
  Fonte secundaria: Fiocruz / Biomanguinhos -- Parceria Tecnologica (Fiocruz PDT)
  Objetivo: eficacia em modelo canino, GLP safety

FASE 3 (30-60 meses) -- Registro e escala: ~R$ 20-50M
  Fonte primaria: Fiocruz / Biomanguinhos -- Parceria Tecnologica (Fiocruz/Biomanguinhos)
  Fonte secundaria: BNDES Saude Animal
  Objetivo: registro MAPA, transferencia de tecnologia

Mercado total enderecavel: ~7,046,000 caes/ano (R$ 880,750,000/ano). O argumento One Health (vacinar caes para proteger humanos) e o diferencial estrategico para justificar investimento publico.

---

*Report generated by the Marley reverse vaccinology pipeline -- Platform A (mRNA optimization)*
