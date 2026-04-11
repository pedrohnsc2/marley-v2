"""09_sae — Sparse Autoencoders para interpretabilidade mecanistica.

Aplica Sparse Autoencoders (SAEs) sobre embeddings de proteinas e RNA
para descobrir features interpretaveis que explicam propriedades
biologicas. Inspirado no trabalho de interpretabilidade mecanistica
de LLMs, adaptado para modelos biologicos.

Intuicao:
    Embeddings do ESM-2 sao densos e nao interpretaveis. SAEs aprendem
    uma representacao esparsa onde cada neuronio ativo corresponde a
    um conceito biologico interpretavel (ex: "dominio transmembrana",
    "motivo de ligacao a zinco", "proteina secretada").

Ref: Cunningham H et al. (2023) Sparse Autoencoders — arXiv:2309.08600
"""

__version__ = "0.1.0"
