"""11_scientist — Motor de sintese cientifica autonomo (AI Scientist).

Integra todos os modulos anteriores (01-10) em um loop de descoberta
deterministico. Cinco agentes especializados leem resultados, cruzam
achados, geram hipoteses priorizadas e propoem experimentos.

Agentes:
    LiteratureAgent   — identifica gaps no conhecimento (RAG + KG)
    DesignAgent       — propoe variantes de ASO/epitopo (EvoDiff + RL)
    ValidationAgent   — avalia candidatos (termodinamica + estrutura)
    KnowledgeAgent    — encontra conexoes inesperadas (KG + ESM + SAE)
    ReportAgent       — sintetiza tudo em narrativa cientifica

Pipeline do loop de descoberta:
    1. Cada agente analisa estado atual dos modulos
    2. Cruzar achados (consenso vs divergencia entre agentes)
    3. Gerar hipoteses priorizadas com evidencias
    4. Propor plano de validacao experimental
    5. Gravar relatorio consolidado

Totalmente deterministico — sem chamadas LLM. Opera como motor
de sintese computacional que cruza dados de 10 modulos.

Ref: Lu C et al. (2024) The AI Scientist — arXiv:2408.06292
"""

__version__ = "0.2.0"
