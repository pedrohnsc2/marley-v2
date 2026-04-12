"""mrl_quantum.vqe — VQE para estrutura eletronica do sitio ativo da RNase H1.

Calcula a energia do estado fundamental do sitio catalítico de dois-metais
da RNase H1 humana (PDB 2QKB) usando Variational Quantum Eigensolver.
O sitio ativo contem dois ions Mg2+ coordenados por residuos D10, E48, D70, H124,
aguas e o substrato RNA:DNA — essencial para a clivagem mediada por ASO.
"""

__all__ = [
    "build_active_site",
    "run_vqe",
    "run_fci_reference",
]
