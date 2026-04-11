"""05_rosettafold — Analise estrutural computacional do complexo ASO:SL RNA.

Modulo de analise 3D que substitui a dependencia de RoseTTAFold2NA
por um pipeline computacional de geometria molecular. Constroi modelos
helicoidais idealizados, analisa parametros estruturais, e decompoe
a energia de ligacao em componentes fisicos.

Funcionalidades:
    - Construcao de duplex RNA:DNA hibrido na forma-A
    - Aplicacao de modificacoes LNA (C3'-endo constrainado)
    - Backbone fosforotioato (PS) com geometria ajustada
    - Analise de sulcos (major/minor groove)
    - Deteccao de pontes de hidrogenio (Watson-Crick + wobble)
    - Estimativa de SASA (Shrake-Rupley simplificado)
    - Avaliacao de acessibilidade da RNase H
    - Decomposicao de energia: stacking, H-bond, eletrostatica, solvatacao
    - Comparacao com PDBs existentes (MRL_ASO_001.pdb, sl_rna.pdb)

O modelo gerado e valido em formato PDB e pode ser visualizado
em PyMOL, ChimeraX, ou Mol* (B-factor codifica: LNA=20, PS=10).

Ref: Saenger W (1984) Principles of Nucleic Acid Structure. Springer.
     Egli M et al. (2005) Chem Biol 12(6):669-675 (LNA geometry).
     Sugimoto N et al. (1995) Biochemistry 34:11211-11216 (RNA:DNA NN params).
"""

__version__ = "1.0.0"
