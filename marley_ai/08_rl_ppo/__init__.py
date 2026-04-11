"""08_rl_ppo — Otimizacao de ASO por Aprendizado por Reforco (REINFORCE).

Usa policy gradient (REINFORCE) para otimizar modificacoes quimicas do
ASO anti-SL RNA de L. infantum. O agente aprende a posicionar LNA,
escolher backbone (PS/PO) e ajustar bases para maximizar eficacia
terapeutica medida por metricas termodinamicas do aso_math.thermo:
    - delta_G de hibridizacao (ligacao ASO:alvo)
    - Tm (estabilidade termica do duplex)
    - GC content (estabilidade biofisica)
    - Hairpin stability (estruturas secundarias indesejaveis)
    - Self-dimer stability (auto-complementaridade)

Implementacao numpy pura — sem PyTorch ou stable-baselines3.

Ref: Williams RJ (1992) Machine Learning 8(3-4):229-256
"""

__version__ = "0.2.0"
