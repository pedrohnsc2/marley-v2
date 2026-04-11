"""04_rna_fm — Analise de RNA com deep learning para o SL RNA de Leishmania.

Combina tres abordagens complementares:
    1. Encoder customizado de RNA (Transformer + Masked Nucleotide Modeling)
    2. Predicao de estrutura secundaria (algoritmo de Nussinov)
    3. Analise de unicidade no espaco de embeddings

Demonstra que o SL RNA de Leishmania e estruturalmente e composicionalmente
distinto dos RNAs do hospedeiro (humano/canino), fundamentando a seletividade
do ASO MRL-ASO-001.

Modelo treinado from scratch em Apple MPS — nao depende de RNA-FM externo.

Refs:
    - Nussinov R et al. (1978) SIAM J Appl Math 35(1):68-82
    - Devlin J et al. (2019) BERT pre-training
    - Chen J et al. (2022) arXiv:2204.00300 (RNA-FM — inspiracao)
"""

__version__ = "1.0.0"
