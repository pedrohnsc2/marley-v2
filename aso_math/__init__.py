"""Suite de validacao matematica para MRL-ASO-001.

Cinco modulos independentes que provam a otimalidade do ASO candidato
contra o Spliced Leader RNA de Leishmania infantum:

1. Paisagem termodinamica (mutantes pontuais + varredura de comprimento)
2. Prova de seletividade (complementaridade maxima + entropia posicional)
3. Conservacao evolutiva (alinhamento cross-species do SL RNA)
4. Otimizacao exaustiva (enumeracao completa + posicoes LNA)
5. Modelo de resistencia (Markov + Poisson para escape por mutacao)

Cada modulo pode ser executado independentemente:
    python -m aso_math.01_thermodynamic_landscape.run
"""

__version__ = "1.0.0"
