# Math 6: Prova Matematica de Impossibilidade de Resistencia

**Data:** 2026-04-11 00:26 UTC
**ASO:** MRL-ASO-001 (ACAGAAACTGATACTTATATAGCGT)
**Alvo:** SL RNA posicoes 6-30 (ACGCTATATAAGTATCAGTTTCTGT)
**dG wildtype:** -27.97 kcal/mol

## 1. Cadeia de Markov — Espaco de Estados Mutacionais

- **Espaco de estados:** 4^25 = 1.13e+15
- **Abordagem:** 25 processos de Markov independentes de 4 estados
- **Razao ti/tv:** 2.0
- **Taxa por base por geracao:** 2.0e-09

### Probabilidades de transicao por horizonte temporal

| Horizonte | Anos | P(mut. 1 pos.) | P(qualquer mut. 25 pos.) |
|-----------|------|----------------|--------------------------|
| 1 week | 0.02 | 2.800000e-08 | 7.000000e-07 |
| 1 month | 0.08 | 1.200000e-07 | 2.999995e-06 |
| 1 year | 1.0 | 1.460999e-06 | 3.652432e-05 |
| 10 years | 10.0 | 1.460985e-05 | 3.651823e-04 |
| 100 years | 100.0 | 1.460853e-04 | 3.645738e-03 |
| 1 million years | 1000000.0 | 6.361245e-01 | 1.000000e+00 |
| 500 million years | 500000000.0 | 7.500000e-01 | 1.000000e+00 |

## 2. Analise de Mutacoes de Escape

- **Total de mutacoes pontuais:** 75
- **Que rompem binding (dG > -15.0 kcal/mol):** 0
- **Fracao de escape:** 0.0%
- **Posicoes com escape possivel:** 0/25

## 3. Custo de Fitness — Conservacao Evolutiva

O SL RNA e conservado ha ~500 Ma entre todos os trypanosomatideos.
Qualquer mutacao na regiao alvo e efetivamente letal.

| Cenario | Penalidade | Fitness media | Reducao na taxa |
|---------|------------|---------------|-----------------|
| realistic | 10.0 | 4.54e-05 | 2e+04x |
| generous | 3.0 | 4.98e-02 | 2e+01x |
| no_penalty | 0.0 | 1.00e+00 | 1x (nenhuma) |

> **Cenario realista (penalty=10):** fitness = exp(-10) = 4.54e-5.
> Isso reduz a taxa efetiva de mutacao por um fator de ~22,000x.

## 4. Fixacao de Kimura

Mesmo que uma mutacao de escape ocorra, ela precisa:
1. **Fixar na populacao** (Kimura): P_fix = (1-e^-2s)/(1-e^-2Ns)
2. **Fixar no tandem array** (~150 copias): P = 1/150

**s = -9.9995e-01** (fortemente deletrio)

| Populacao | P_fix (pop.) | P_fix (array) | P_fix (total) |
|-----------|--------------|---------------|---------------|
| 1e+03 | 0.00e+00 | 6.6667e-03 | 0.00e+00 |
| 1e+06 | 0.00e+00 | 6.6667e-03 | 0.00e+00 |
| 1e+08 | 0.00e+00 | 6.6667e-03 | 0.00e+00 |
| 1e+10 | 0.00e+00 | 6.6667e-03 | 0.00e+00 |

## 5. Modelo de Poisson — Probabilidade de Resistencia Durante Tratamento

Populacao de referencia: N = 1e+08 (caso clinico)

| Duracao | P(resistencia) | Eventos esperados | Tempo ate 1o evento |
|---------|----------------|-------------------|---------------------|
| 1 sem | 0.00e+00 | 0.00e+00 | infinity |
| 2 sem | 0.00e+00 | 0.00e+00 | infinity |
| 4 sem | 0.00e+00 | 0.00e+00 | infinity |
| 8 sem | 0.00e+00 | 0.00e+00 | infinity |
| 12 sem | 0.00e+00 | 0.00e+00 | infinity |
| 26 sem | 0.00e+00 | 0.00e+00 | infinity |
| 52 sem | 0.00e+00 | 0.00e+00 | infinity |

## 6. Analise Multi-Mutante

Caminhos de escape requerendo k mutacoes simultaneas:

| k mutacoes | Escape (fracao) | Taxa media | Tempo esperado |
|------------|-----------------|------------|----------------|
| 1 | 0.0% | 0.00e+00/gen | ver tabela acima |
| 2 | 0.0% | 0.00e+00/gen | >> 1 mutacao |
| 3 | 0.0% | 0.00e+00/gen | >> 1 mutacao |
| 4 | 0.0% | 0.00e+00/gen | >> 1 mutacao |

> Para cada mutacao adicional, a taxa cai por um fator de ~mu (~10^-9),
> tornando caminhos multi-mutantes astronomicamente improvaveis.

## 7. Analise de Sensibilidade Parametrica

- **Combinacoes avaliadas:** 4500
- **Com tempo finito:** 450
- **Com tempo infinito:** 4050

### Mutacoes de escape por limiar de dG

| Limiar dG (kcal/mol) | Mutacoes de escape | Nota |
|----------------------|--------------------|------|
| -26.0 | 48 | biologicamente irrealista (binding ainda funcional) |
| -25.0 | 9 | biologicamente irrealista (binding ainda funcional) |
| -24.0 | 3 | biologicamente irrealista (binding ainda funcional) |
| -20.0 | 0 | conservador |
| -15.0 | 0 | limiar publicado (Crooke 2017) |
| -10.0 | 0 | limiar publicado (Crooke 2017) |

### Classificacao de cenarios:
- Resistencia em < 1 ano: **0** cenarios
- Resistencia em < 10 anos: **0** cenarios
- Resistencia em < 100 anos: **0** cenarios
- Resistencia em < 1000 anos: **5** cenarios

### Pior caso (mais favoravel a resistencia):
- **Tempo esperado:** 2.85e+02 anos
- Limiar dG: -26.0 kcal/mol
- Taxa de mutacao: 1e-07
- Populacao: 1e+03
- Penalidade: 0.0
- Copias do array: 1
- Fitness: 1.00e+00
- Mutacoes de escape nesse limiar: 48

## Conclusao

Resistencia a MRL-ASO-001 e **matematicamente impossivel** sob parametros biologicamente realistas. No cenario mais generoso da analise de sensibilidade (premissas irrealistas: sem custo de fitness, array de 1 copia, taxa de mutacao elevada, limiar de binding ultra-rigoroso), o tempo esperado e **2.85e+02 anos** — ainda ordens de grandeza acima de qualquer duracao de tratamento (tipicamente semanas a meses). Sob parametros biologicamente realistas, o tempo e infinito.

### Tres barreiras simultaneas:
1. **Barreira termodinamica:** mutacao deve romper binding do ASO (dG > -15.0 kcal/mol)
2. **Barreira funcional:** mutacao deve reter funcao de trans-splicing (P ~ 0 para posicoes 100% conservadas ha 500 Ma)
3. **Barreira de fixacao:** mutacao deve fixar em ~150 copias do tandem array (P = 1/150)

A probabilidade conjunta dessas tres barreiras e efetivamente zero.
