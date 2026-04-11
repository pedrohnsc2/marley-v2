"""02_leish_kg — Knowledge Graph de Leishmania.

Constroi um grafo de conhecimento integrando dados de:
    - TriTrypDB (genes, proteinas, vias metabolicas)
    - UniProt (anotacoes funcionais, GO terms)
    - PubMed (relacoes extraidas de abstracts via NLP)
    - Resultados do pipeline Marley (epitopos, alvos ASO)

O grafo conecta entidades biologicas (proteinas, genes, vias,
farmacos, epitopos) por relacoes tipadas, permitindo raciocinio
sobre o parasita e descoberta de novos alvos terapeuticos.
"""

__version__ = "0.1.0"
