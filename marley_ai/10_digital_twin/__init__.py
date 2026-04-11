"""10_digital_twin — Digital twin imunologico canino.

Simula a resposta imune de um cao a construtos vacinais e ASOs
terapeuticos. Modela compartimentos biologicos (injecao, linfonodo,
baco, figado) e celulas imunes (APCs, T CD4+, T CD8+, B cells)
como um sistema de equacoes diferenciais informado por dados.

O digital twin permite:
    - Predizer resposta imune antes de experimentos in vivo
    - Otimizar protocolo de vacinacao (dose, intervalos, adjuvantes)
    - Simular cenarios de co-infeccao (Leishmania + Erlichia, etc.)
    - Estimar janela terapeutica para ASO

Ref: Laubenbacher R et al. (2022) Digital twins in medicine — Nat Comp Sci
"""

__version__ = "0.1.0"
