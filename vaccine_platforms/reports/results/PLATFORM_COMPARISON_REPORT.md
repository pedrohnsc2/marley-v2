# Cross-Platform Comparison Report -- Marley Vaccine Pipeline

**Generated:** 2026-04-11 01:39 UTC
**Pipeline version:** 0.1.0
**Platforms compared:** 3 (mRNA-LNP, E. coli recombinante, L. tarentolae LEXSY)
**Target:** Leishmania infantum -- Leishmaniose Visceral Canina (LVC)

---

## 1. Comparison Matrix

### 1.1 Construct Properties

| Metric | Platform A (mRNA) | Platform B (E. coli) | Platform C (L. tarentolae) |
|---|---|---|---|
| Comprimento da proteina (aa) | 287 | 278 | 302 |
| Peso molecular (kDa) | 30.7 | 30.2 | 32.6 |
| Indice de instabilidade | 32.02 | 38.02 | 30.58 |
| Proteina estavel | Sim | Sim | Sim |
| Epitopos unicos | 11 | 11 | 11 |
| Proxy de antigenicidade (0-1) | 0.6245 | 0.6245 | 0.6557 |

### 1.2 Production & Expression

| Metric | Platform A (mRNA) | Platform B (E. coli) | Platform C (L. tarentolae) |
|---|---|---|---|
| Sistema de expressao | In vivo (celulas do hospedeiro) | E. coli BL21(DE3) / pET-28a(+) | L. tarentolae LEXSY / pLEXSY-sat2 |
| Purificacao | N/A (mRNA encapsulado em LNP) | IMAC (Ni-NTA) + TEV + IMAC reverso + SEC | Secretada + His6/Strep-tag tandem (ou vacina viva sem purificacao) |
| Adjuvante | LNP (auto-adjuvante, SM-102:DSPC:Chol:PEG-DMG) | QuilA saponina ($1.50/dose) | Intrinseco (LPG, GP63, GIPLs -- TLR2/TLR4) |
| Modificacoes pos-traducionais | Nao (traducao in vivo) | Nao (procariotico) | Sim (N-glicosilacao, GPI, acetilacao N-terminal) |

### 1.3 Logistics & Infrastructure

| Metric | Platform A (mRNA) | Platform B (E. coli) | Platform C (L. tarentolae) |
|---|---|---|---|
| Cadeia de frio | -70C (ou liofilizado 2-8C) | 2-8C | 2-8C (viva) ou temp. ambiente (liofilizada) |
| Infraestrutura necessaria | Alta (GMP mRNA, LNP, NTPs modificados) | Baixa (shaker flask, BL21(DE3)) | Baixa (BSL-1, 26C, sem CO2) |

### 1.4 Cost per Dose

| Scale | Platform A (mRNA) | Platform B (E. coli) | Platform C Live | Platform C Protein |
|---|---|---|---|---|
| Lab | $260.01 | $4.79 | $0.00066 | $3.07 |
| Industrial | $5.75 | $2.28 | $0.00008 | $2.13 |

### 1.5 ASO Compatibility

| Platform | Compatible | Details |
|---|---|---|
| A (mRNA) | YES | Totalmente compativel -- mRNA nao depende de SL RNA |
| B (E. coli) | YES | Totalmente compativel -- E. coli nao depende de SL RNA |
| C (L. tarentolae) | **WARNING -- NO** | WARNING -- DeltaDeltaG = 0.00 kcal/mol (limiar: 2.0). O ASO MRL-ASO-001 tem complementaridade perfeita com o SL RNA de L. tarentolae. Requer separacao temporal ou modificacao do ASO. |

---

## 2. Radar Chart Data (6-Axis, 0-1 Normalized)

```
Axis                    | Platform A (mRNA) | Platform B (E. coli) | Platform C (L. tarentolae)
Eficacia (antigenicidade)| 0.6245            | 0.6245               | 0.6557
Eficiencia de custo     | 0.0               | 0.6037               | 1.0
Viabilidade logistica   | 0.3               | 0.85                 | 0.9
Facilidade de manufatura| 0.35              | 0.8                  | 0.9
Prontidao regulatoria   | 0.3               | 0.75                 | 0.55
Compatibilidade ASO     | 1.0               | 1.0                  | 0.0
```

**Scoring methodology:**
- **efficacy:** Proxy de antigenicidade do construto (0-1). L. tarentolae recebe bonus de 5% por glicosilacao nativa.
- **cost_efficiency:** 1 - (custo_industrial / max_custo). Custo menor = score maior.
- **logistic_viability:** Media de cold_chain_score e infrastructure_score. Cold chain: -70C=0.3, -20C=0.5, 2-8C=0.8, temp.amb=1.0. Infraestrutura: alta=0.3, media=0.6, baixa=0.9.
- **manufacturing_ease:** Score qualitativo: mRNA (IVT+LNP)=0.35, E.coli (fermentacao+IMAC)=0.80, L.tarentolae viva (cultura)=0.90.
- **regulatory_readiness:** Score qualitativo baseado em precedentes regulatorios: mRNA vet=0.30 (sem precedente), E.coli recomb=0.75 (Leish-Tec como precedente), L.tarentolae viva=0.55 (vetor vivo, sem precedente especifico).
- **aso_compatibility:** Binario: 1.0 se compativel com MRL-ASO-001, 0.0 se incompativel (issue de ortogonalidade SL RNA).

---

## 3. Cost Comparison

### 3.1 Per Dose (USD)

| Scale | Platform A | Platform B | Platform C (live) | Platform C (protein) | Leish-Tec |
|---|---|---|---|---|---|
| Lab | $260.01 | $4.79 | $0.00066 | $3.07 | $27.03 |
| Industrial | $5.75 | $2.28 | $0.00008 | $2.13 | $27.03 |

### 3.2 Per Treatment Course (USD)

*Plataforma A: 2 doses; Plataformas B e C: 3 doses; Leish-Tec: 3 doses + reforco anual*

| Scale | Platform A (2 doses) | Platform B (3 doses) | Platform C live (3 doses) | Platform C protein (3 doses) | Leish-Tec (3 doses) |
|---|---|---|---|---|---|
| Lab | $520.02 | $14.36 | $0.00198 | $9.20 | $81.09 |
| Industrial | $11.50 | $6.84 | $0.00024 | $6.38 | $81.09 |

### 3.3 Cost Reduction vs Leish-Tec (Industrial)

| Platform | Reduction (%) |
|---|---|
| A (mRNA) | 85.8% |
| B (E. coli) | 91.6% |
| C (live) | 100.0% |
| C (protein) | 92.1% |

---

## 4. Timeline Estimates

| Milestone | Platform A (mRNA) | Platform B (E. coli) | Platform C (L. tarentolae) |
|---|---|---|---|
| First testable batch | 6 months | 3 months | 4 months |
| Pre-clinical data | 18 months | 12 months | 15 months |
| Regulatory submission | +12 months | +9 months | +12 months |
| **Total to market** | **36-60 meses** | **24-36 meses** | **31-48 meses** |
| Bottleneck | Via regulatoria: primeiro mRNA veterinario no Brasil. Requer... | Competicao direta com Leish-Tec. Diferencial: multi-epitopo ... | Ortogonalidade com ASO limita uso combinado. Regulatorio par... |

---

## 5. Strategic Recommendation

> As tres plataformas nao sao mutuamente exclusivas. Cada uma atende a um segmento diferente de mercado, financiamento e estrategia regulatoria. A recomendacao e prosseguir com as tres em paralelo, priorizando Platform B para entrada rapida no mercado veterinario e Platform A para soberania vacinal.

### 5.1 Platform Roles

#### Platform A -- mRNA-LNP
- **Strategic role:** Soberania vacinal + financiamento governamental
- **Target market:** Programa Nacional de Controle da LV (SUS/MAPA), campanhas de vacinacao em massa
- **Funding:** Fiocruz/Bio-Manguinhos (planta mRNA), FINEP Inovacao, BNDES Saude Animal
- **Key advantage:** Plataforma de mRNA ja sendo instalada no Brasil (Fiocruz + BioNTech). Demonstrar dominio da tecnologia para alem de COVID-19.
- **Key risk:** Via regulatoria incerta, custo mais alto, cold chain exigente (-70C ou liofilizacao)
- **Priority:** ALTA -- argumento de soberania vacinal e forte para financiamento publico
- **ASO compatibility:** Totalmente compativel

#### Platform B -- E. coli Recombinante
- **Strategic role:** Entrada rapida no mercado veterinario
- **Target market:** Clinicas veterinarias privadas, mercado pet premium
- **Funding:** Investimento privado, PIPE-FAPESP, parceria com industria veterinaria (Ourofino, Ceva)
- **Key advantage:** Menor custo ($2.28/dose industrial), infraestrutura simples, via regulatoria conhecida (precedente Leish-Tec). Prazo mais curto ate o mercado (24-36 meses).
- **Key risk:** Competicao direta com Leish-Tec. Sem modificacoes pos-traducionais (sem glicosilacao).
- **Priority:** ALTA -- via mais rapida e barata para comercializacao
- **ASO compatibility:** Totalmente compativel

#### Platform C -- L. tarentolae LEXSY
- **Strategic role:** Parceria academica + diferencial cientifico
- **Target market:** Colaboracao Fiocruz, mercado de alta performance (glicosilacao nativa)
- **Funding:** MCTI/CNPq, parceria academica com grupos de Leishmania, eventual Fiocruz PDT
- **Key advantage:** Unica plataforma com modificacoes pos-traducionais nativas de Leishmania (glicosilacao, GPI). Adjuvante intrinseco elimina custo de adjuvante exogeno. Custo minimo ($0.00008/dose viva, industrial).
- **Key risk:** CRITICO: Problema de ortogonalidade com ASO MRL-ASO-001. O ASO tem complementaridade perfeita com o SL RNA de L. tarentolae (DeltaDeltaG = 0.00 kcal/mol). Uso combinado ASO + vacina C requer separacao temporal rigorosa ou modificacao do ASO.
- **Priority:** MEDIA -- diferencial cientifico forte, mas issue de ortogonalidade limita uso combinado
- **ASO compatibility:** INCOMPATIVEL sem precaucoes (requer separacao temporal ou modificacao do ASO)

**Mitigation options for ASO orthogonality issue:**
1. Administrar vacina L. tarentolae ANTES do tratamento com ASO (intervalo de 2-4 semanas)
1. Usar intervalo longo (>1 semana) entre ASO e vacinacao
1. Modificar ASO com LNA/2'-OMe para aumentar seletividade por L. infantum
1. Considerar Platform B (E. coli) quando compatibilidade ASO e prioridade

### 5.2 Combined Strategy (Three Phases)

| Phase | Timeline | Budget | Primary | Secondary | Funding |
|---|---|---|---|---|---|
| Fase 1: Prova de Conceito (0-12 meses) | 0-12 mo | R$ 500K - 2M | B (E. coli) -- dados pre-clinicos em camundongos | A (mRNA) -- prova de conceito em modelo murino | MCTI/CNPq + FAP estadual |
| Fase 2: Pre-clinico Avancado (12-30 meses) | 12-30 mo | R$ 5M - 10M | B (E. coli) -- estudo em caes, comparacao com Leish-Tec | A (mRNA) -- estudo de imunogenicidade em caes | FINEP Inovacao + Fiocruz PDT |
| Fase 3: Registro e Escala (30-60 meses) | 30-60 mo | R$ 20M - 50M | B (E. coli) -- registro MAPA, producao piloto | A (mRNA) -- submissao regulatoria (via de novo) | Fiocruz/Bio-Manguinhos + BNDES |

### 5.3 Market Context

| Metric | Value |
|---|---|
| Total dogs in Brazil | 54,200,000 |
| Dogs at risk (endemic areas) | 35,230,000 |
| Addressable market (doses/year) | 14,092,000 |
| Market value (BRL/year) | R$ 880,750,000 |
| Competitor | Leish-Tec (Hertape-Calier, proteina A2 + saponina, ~R$150/dose) |
| Marley differential | Multi-epitopo (5 alvos proteicos vs 1 alvo), pipeline computacional completo com ASO terapeutico complementar |

---

*Report generated by the Marley reverse vaccinology pipeline -- Cross-Platform Comparison Module*
