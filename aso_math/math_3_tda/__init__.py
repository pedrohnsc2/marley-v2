"""Math 3 — Analise de Dados Topologicos (TDA) de estruturas ASO e SL RNA.

Aplica Topological Data Analysis para caracterizar a topologia das
estruturas 3D do ASO livre, SL RNA alvo e duplex ASO:SL RNA.

Implementacao:
- Extracao de nuvens de pontos de arquivos PDB (backbone atoms)
- Filtracao de Vietoris-Rips com union-find para H0
- Diagramas de persistencia para H0 e H1
- Distancia bottleneck entre diagramas
- Score de estabilidade topologica
"""
