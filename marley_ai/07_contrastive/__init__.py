"""07_contrastive — Aprendizado contrastivo epitopo-MHC canino.

Treina um modelo de aprendizado contrastivo (estilo CLIP) para alinhar
representacoes de epitopos e alelos DLA no mesmo espaco vetorial.
Permite predizer afinidade de ligacao para pares epitopo-DLA nao vistos
e identificar novos epitopos candidatos.

Intuicao biologica:
    Epitopos que ligam ao mesmo alelo DLA devem ficar proximos no espaco
    de embeddings; epitopos que nao ligam devem ficar distantes.
    Isso aprende uma "funcao de afinidade" implicita sem precisar de
    dados estruturais de cristalografia.

Ref: Radford A et al. (2021) CLIP — ICML 2021
"""

__version__ = "0.1.0"
