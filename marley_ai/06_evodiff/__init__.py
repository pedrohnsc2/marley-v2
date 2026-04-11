"""06_evodiff — Geracao de sequencias via difusao discreta.

Implementa um modelo de difusao discreta (D3PM simplificado) treinado
do zero em PyTorch para gerar variantes de ASOs e epitopos com
propriedades otimizadas para o projeto Marley.

Dois modos de operacao:
    - ASO: gera variantes de nucleotideos a partir de MRL-ASO-001
    - Epitopo: gera peptideos candidatos a partir dos 11 epitopos canonicos

O processo de difusao corrompe sequencias gradualmente (forward) e
aprende a reconstrui-las (reverse), permitindo gerar novas sequencias
a partir de ruido puro.

Ref: Alamdari S et al. (2023) EvoDiff — NeurIPS 2023
     Austin J et al. (2021) Structured Denoising Diffusion — NeurIPS 2021
"""

__version__ = "0.2.0"
