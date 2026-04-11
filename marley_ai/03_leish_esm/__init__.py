"""03_leish_esm — Embeddings proteicos via ESM-2 para Leishmania.

Gera representacoes vetoriais densas do proteoma de L. infantum usando
o modelo ESM-2 (Evolutionary Scale Modeling). Esses embeddings capturam
propriedades estruturais e funcionais aprendidas de milhoes de proteinas,
permitindo:
    - Busca de homologos funcionais sem BLAST
    - Clustering de proteinas por funcao predita
    - Features para o modelo contrastivo (modulo 07)
    - Input para o digital twin (modulo 10)

Ref: Lin Z et al. (2023) Science 379(6637):1123-1130
"""

__version__ = "0.1.0"
