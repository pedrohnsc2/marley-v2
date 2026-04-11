"""Pipeline de entrega e farmacocinetica para MRL-ASO-001.

Seis modulos que modelam o caminho do ASO desde a formulacao ate
a acao terapeutica no macrofago infectado por Leishmania:

A. Estabilidade quimica (degradacao em soro, meia-vida, modificacoes)
B. Permeabilidade de membrana (uptake celular, endosomal escape)
C. Conjugados de entrega (GalNAc, peptideos, anticorpos)
D. Nanoparticulas lipidicas (formulacao LNP, parametros criticos)
E. ADMET in silico (absorcao, distribuicao, metabolismo, excrecao, toxicidade)
F. Modelo imune estocastico (SDE para resposta imune inata ao ASO)

Cada modulo pode ser executado independentemente:
    python -m aso_delivery.module_a_stability.run
"""

__version__ = "0.1.0"
