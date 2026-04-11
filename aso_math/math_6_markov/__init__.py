"""Math 6 — Prova matematica de impossibilidade de resistencia via cadeia de Markov.

Combina modelo de Markov para transicoes mutacionais, custos de fitness
baseados em conservacao, modelo de Poisson para emergencia de resistencia,
e probabilidade de fixacao de Kimura para demonstrar que resistencia a
MRL-ASO-001 e efetivamente impossivel sob qualquer cenario biologicamente
plausivel.

Estende o modelo existente em 05_resistance_model com:
- Cadeia de Markov para 25 posicoes independentes (transicoes vs transversoes)
- Custos de fitness derivados de conservacao evolutiva (~500 Ma)
- Modelo de Poisson para emergencia em populacao finita
- Fixacao de Kimura com selecao negativa
- Analise multi-mutante para caminhos de escape de 2+ mutacoes
"""
