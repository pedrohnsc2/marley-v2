"""A4 -- Analise estrategica de financiamento e mercado.

Gera uma analise estruturada do contexto estrategico para a vacina
de mRNA contra leishmaniose visceral canina no Brasil:

    - Programas de fomento do governo brasileiro
    - Argumento de soberania vacinal (precedente COVID-19)
    - Parceria com Fiocruz/Biomanguinhos
    - Via regulatoria USDA para biologicos veterinarios
    - Tamanho de mercado e populacao canina em risco
    - Scoring de cada fonte de financiamento

Os dados de mercado sao baseados em fontes oficiais (IBGE, SINAN,
MAPA) e literatura recente sobre leishmaniose visceral no Brasil.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

# ---------------------------------------------------------------------------
# Dados de mercado -- Brasil
# ---------------------------------------------------------------------------

# Populacao canina no Brasil (IBGE, Pesquisa Nacional de Saude 2019)
# ~54.2 milhoes de caes domiciliados
TOTAL_DOGS_BRAZIL: Final[int] = 54_200_000

# Caes em areas endemicas de LV (estimativa conservadora)
# Regioes Norte, Nordeste e partes do Sudeste/Centro-Oeste
# ~65% da populacao canina esta em areas com transmissao ativa
DOGS_AT_RISK_FRACTION: Final[float] = 0.65
DOGS_AT_RISK: Final[int] = round(TOTAL_DOGS_BRAZIL * DOGS_AT_RISK_FRACTION)

# Prevalencia media de LVC em areas endemicas (soro-prevalencia canina)
# Varia de 3% a 75% dependendo da regiao (meta-analise)
SEROPREVALENCE_RANGE: Final[tuple[float, float]] = (0.03, 0.75)
SEROPREVALENCE_MEAN: Final[float] = 0.25

# Casos humanos notificados por ano (SINAN/MS, media 2018-2022)
HUMAN_CASES_PER_YEAR: Final[int] = 3_500
HUMAN_CASE_FATALITY_RATE: Final[float] = 0.08  # ~8% letalidade

# Custo do programa de controle de LV no Brasil (estimativa)
# Inclui: diagnostico canino, eutanasia de soropositivos, borrifacao,
# tratamento humano, vigilancia epidemiologica
ANNUAL_CONTROL_COST_BRL: Final[float] = 350_000_000.0  # R$ 350 milhoes/ano

# Preco de vacinas veterinarias no Brasil (referencia)
LEISH_TEC_PRICE_BRL: Final[float] = 250.0  # R$ 250 por 3 doses (Hertape-Calier)

# Cambio referencia (abril 2026)
BRL_TO_USD: Final[float] = 0.18  # 1 BRL ~ 0.18 USD


# ---------------------------------------------------------------------------
# Fontes de financiamento
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FundingSource:
    """Fonte de financiamento potencial para o projeto."""

    name: str
    institution: str
    description: str
    typical_amount_brl: tuple[float, float]  # (min, max) em BRL
    probability_pct: float  # 0-100
    timeline_months: tuple[int, int]  # (min, max) meses ate o recurso
    requirements: list[str]
    strategic_fit: str
    score: float  # 0-10, pontuacao geral


# Fontes de financiamento avaliadas
FUNDING_SOURCES: Final[tuple[FundingSource, ...]] = (
    FundingSource(
        name="BNDES Saude Animal",
        institution="Banco Nacional de Desenvolvimento Economico e Social",
        description=(
            "Linha de credito para inovacao em saude animal. "
            "Prioridade para projetos que reduzam dependencia de importacoes "
            "e desenvolvam tecnologia nacional."
        ),
        typical_amount_brl=(2_000_000.0, 20_000_000.0),
        probability_pct=35.0,
        timeline_months=(12, 24),
        requirements=[
            "Empresa brasileira como proponente (parceria academia-industria)",
            "Plano de negocio com viabilidade economica",
            "Contrapartida financeira da empresa (20-40%)",
            "Producao nacional obrigatoria",
        ],
        strategic_fit=(
            "Forte alinhamento: vacina veterinaria com tecnologia nacional "
            "de mRNA, reducao de importacao de biologicos. Contudo, exige "
            "parceiro industrial ja estabelecido."
        ),
        score=6.5,
    ),
    FundingSource(
        name="FINEP Inovacao",
        institution="Financiadora de Estudos e Projetos",
        description=(
            "Subvencao economica e credito para inovacao tecnologica. "
            "Programa Tecnova e Chamadas Publicas de Subvencao Economica "
            "financiam P&D em saude animal e biotecnologia."
        ),
        typical_amount_brl=(500_000.0, 10_000_000.0),
        probability_pct=40.0,
        timeline_months=(8, 18),
        requirements=[
            "Projeto de inovacao tecnologica com diferencial claro",
            "Equipe tecnica qualificada (doutores, patentes)",
            "Infraestrutura laboratorial comprovada",
            "Potencial de mercado documentado",
        ],
        strategic_fit=(
            "Excelente: FINEP historicamente financia vacinas e biotecnologia. "
            "A tecnologia de mRNA se enquadra em 'inovacao de fronteira'. "
            "Subvencao economica (nao-reembolsavel) e ideal para fase pre-clinica."
        ),
        score=7.5,
    ),
    FundingSource(
        name="MCTI / CNPq -- Chamada Universal + INCT",
        institution="Ministerio da Ciencia, Tecnologia e Inovacao",
        description=(
            "Financiamento academico via editais universais do CNPq e "
            "Institutos Nacionais de Ciencia e Tecnologia (INCTs). "
            "Valores menores, mas sem exigencia de parceiro industrial."
        ),
        typical_amount_brl=(100_000.0, 2_000_000.0),
        probability_pct=55.0,
        timeline_months=(6, 12),
        requirements=[
            "Pesquisador com vinculo institucional (universidade/instituto)",
            "Curriculo Lattes atualizado com producao relevante",
            "Projeto com impacto cientifico e social",
            "Parceria inter-institucional (bonus na avaliacao)",
        ],
        strategic_fit=(
            "Muito bom para fase de prova de conceito e estudos pre-clinicos. "
            "Valores menores limitam escala, mas a alta probabilidade de "
            "aprovacao e a flexibilidade de uso compensam."
        ),
        score=7.0,
    ),
    FundingSource(
        name="FAPESP / FAPs estaduais -- PIPE/PAPPE",
        institution="Fundacoes Estaduais de Amparo a Pesquisa",
        description=(
            "Programas de inovacao para pequenas empresas (PIPE/SP, "
            "PAPPE/outros estados). Financiam P&D em parceria com ICTs. "
            "Especialmente forte em SP (FAPESP PIPE) e MG (FAPEMIG)."
        ),
        typical_amount_brl=(200_000.0, 2_500_000.0),
        probability_pct=45.0,
        timeline_months=(6, 18),
        requirements=[
            "Empresa de pequeno porte (EPP) como proponente",
            "Sede no estado da FAP financiadora",
            "Projeto com potencial de produto comercializavel",
            "Pesquisador vinculado a ICT como co-PI",
        ],
        strategic_fit=(
            "Bom para MVP e prova de conceito. PIPE-FAPESP tem historico "
            "de apoiar startups de biotecnologia. Ideal se houver spin-off "
            "academico ou empresa incubada."
        ),
        score=6.0,
    ),
    FundingSource(
        name="Fiocruz / Biomanguinhos -- Parceria Tecnologica",
        institution="Fundacao Oswaldo Cruz / Bio-Manguinhos",
        description=(
            "Parceria de desenvolvimento tecnologico (PDT) com Bio-Manguinhos. "
            "Biomanguinhos e o maior produtor de vacinas do Brasil e esta "
            "investindo em plataforma de mRNA desde o COVID-19."
        ),
        typical_amount_brl=(5_000_000.0, 50_000_000.0),
        probability_pct=20.0,
        timeline_months=(18, 36),
        requirements=[
            "Tecnologia com potencial de transferencia para producao",
            "Dados pre-clinicos robustos (eficacia em modelo animal)",
            "Alinhamento com portfolio de Biomanguinhos",
            "Demanda de saude publica comprovada",
        ],
        strategic_fit=(
            "O melhor cenario de longo prazo: Biomanguinhos instalou planta "
            "de mRNA em 2023 (acordo com BioNTech) para produzir vacinas "
            "de mRNA no Brasil. A leishmaniose e uma das prioridades do "
            "SUS. Contudo, o foco historico e saude humana -- seria necessario "
            "argumentar o controle da LVC como estrategia One Health para "
            "reduzir casos humanos."
        ),
        score=8.0,
    ),
    FundingSource(
        name="USDA SBIR / International Collaborations",
        institution="US Department of Agriculture",
        description=(
            "Small Business Innovation Research (SBIR) grants para "
            "biologicos veterinarios inovadores. Aberto a colaboracoes "
            "internacionais em certos programas."
        ),
        typical_amount_brl=(750_000.0, 5_000_000.0),
        probability_pct=15.0,
        timeline_months=(12, 24),
        requirements=[
            "Empresa ou instituicao com presenca nos EUA (ou colaborador)",
            "Proposta em ingles com dados preliminares",
            "Via regulatoria USDA-CVB mapeada",
            "Potencial de comercializacao nos EUA",
        ],
        strategic_fit=(
            "Nicho: LV nao e endemica nos EUA, mas importacao de caes "
            "infectados do Mediterraneo e da America Latina e preocupacao "
            "crescente. O USDA-APHIS emitiu alertas sobre L. infantum em "
            "foxhounds americanos. Argumento de biosseguranca valido, "
            "mas a probabilidade e mais baixa."
        ),
        score=4.5,
    ),
)


# ---------------------------------------------------------------------------
# Dataclasses de resultado
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MarketAnalysis:
    """Analise de tamanho de mercado."""

    total_dogs_brazil: int
    dogs_at_risk: int
    dogs_at_risk_fraction: float
    seroprevalence_mean: float
    estimated_infected_dogs: int
    addressable_market_dogs: int
    addressable_market_doses: int
    market_value_brl_year: float
    market_value_usd_year: float
    human_cases_per_year: int
    human_fatality_rate: float
    annual_control_cost_brl: float
    competitor_price_brl: float


@dataclass(frozen=True, slots=True)
class SovereigntyArgument:
    """Argumento de soberania vacinal para mRNA."""

    title: str
    context: str
    precedent: str
    brazil_mrna_capacity: str
    one_health_argument: str
    strategic_value: str


@dataclass
class StrategicAnalysis:
    """Resultado completo da analise estrategica."""

    market: MarketAnalysis
    sovereignty: SovereigntyArgument
    funding_sources: list[dict]
    funding_ranked: list[str]
    regulatory_usda: dict
    regulatory_brazil: dict
    recommendation: str


# ---------------------------------------------------------------------------
# Funcoes de calculo
# ---------------------------------------------------------------------------


def _analyze_market() -> MarketAnalysis:
    """Analisa o mercado potencial para a vacina de mRNA contra LVC.

    Estima:
    - Populacao canina em risco (~35 milhoes)
    - Mercado enderecavel (caes em areas endemicas, proprietarios dispostos a vacinar)
    - Valor de mercado anual

    Returns:
        MarketAnalysis com dados de mercado
    """
    estimated_infected = round(DOGS_AT_RISK * SEROPREVALENCE_MEAN)

    # Mercado enderecavel: 20% dos caes em risco (proprietarios com
    # capacidade e disposicao de vacinar, acesso a veterinario)
    addressable_fraction = 0.20
    addressable_dogs = round(DOGS_AT_RISK * addressable_fraction)

    # Doses por ano (regime de 2 doses/animal, campanha anual)
    doses_per_year = addressable_dogs * 2

    # Valor de mercado (preco competitivo: metade da Leish-Tec)
    target_price_brl = LEISH_TEC_PRICE_BRL * 0.5  # R$ 125 por regime
    market_value_brl = addressable_dogs * target_price_brl
    market_value_usd = market_value_brl * BRL_TO_USD

    return MarketAnalysis(
        total_dogs_brazil=TOTAL_DOGS_BRAZIL,
        dogs_at_risk=DOGS_AT_RISK,
        dogs_at_risk_fraction=DOGS_AT_RISK_FRACTION,
        seroprevalence_mean=SEROPREVALENCE_MEAN,
        estimated_infected_dogs=estimated_infected,
        addressable_market_dogs=addressable_dogs,
        addressable_market_doses=doses_per_year,
        market_value_brl_year=round(market_value_brl, 2),
        market_value_usd_year=round(market_value_usd, 2),
        human_cases_per_year=HUMAN_CASES_PER_YEAR,
        human_fatality_rate=HUMAN_CASE_FATALITY_RATE,
        annual_control_cost_brl=ANNUAL_CONTROL_COST_BRL,
        competitor_price_brl=LEISH_TEC_PRICE_BRL,
    )


def _sovereignty_argument() -> SovereigntyArgument:
    """Constroi o argumento de soberania vacinal para mRNA veterinario.

    O precedente COVID-19 demonstrou a importancia de capacidade
    nacional de producao de vacinas de mRNA. O Brasil investiu
    pesadamente em infraestrutura pos-pandemia (planta de mRNA em
    Bio-Manguinhos). Estender para veterinaria e o proximo passo.

    Returns:
        SovereigntyArgument estruturado
    """
    return SovereigntyArgument(
        title="Soberania vacinal em mRNA: do COVID-19 a One Health",
        context=(
            "A pandemia de COVID-19 expôs a dependencia do Brasil de "
            "tecnologia estrangeira para vacinas de mRNA. O pais foi um "
            "dos ultimos a acessar vacinas de mRNA, com impacto direto em "
            "mortalidade. Desde entao, o governo brasileiro investiu em "
            "capacidade nacional de mRNA."
        ),
        precedent=(
            "Em 2023, Bio-Manguinhos/Fiocruz firmou acordo de transferencia "
            "de tecnologia com BioNTech para producao de vacinas de mRNA no "
            "Brasil. A planta de producao em Santa Cruz (RJ) esta em fase "
            "de instalacao. O investimento total excede R$ 1 bilhao."
        ),
        brazil_mrna_capacity=(
            "A infraestrutura de mRNA em instalacao no Brasil pode ser "
            "utilizada para multiplos produtos alem de COVID-19. Vacinas "
            "veterinarias de mRNA sao candidatas naturais: mesma plataforma "
            "tecnologica, menor barreira regulatoria que saude humana, "
            "demanda interna massiva (leishmaniose e endemica no pais)."
        ),
        one_health_argument=(
            "A leishmaniose visceral e uma zoonose: o cao e o principal "
            "reservatorio urbano de Leishmania infantum. Vacinar caes "
            "reduz a prevalencia canina e, consequentemente, a transmissao "
            "para humanos. O PNCLV (Programa Nacional de Controle da LV) "
            "do Ministerio da Saude atualmente depende de eutanasia de caes "
            "soropositivos -- estrategia eticamente controversa e de "
            "eficacia limitada. Uma vacina eficaz transformaria o programa."
        ),
        strategic_value=(
            "Desenvolver a primeira vacina de mRNA veterinaria do Brasil: "
            "(1) demonstra dominio da plataforma de mRNA para alem de COVID, "
            "(2) cria precedente regulatorio para mRNA veterinario no MAPA, "
            "(3) gera propriedade intelectual nacional, "
            "(4) reduz gastos publicos com controle de LV (~R$ 350 M/ano), "
            "(5) posiciona o Brasil como lider em One Health no G20."
        ),
    )


def _regulatory_usda() -> dict:
    """Via regulatoria USDA para biologicos veterinarios.

    O USDA-APHIS Center for Veterinary Biologics (CVB) regulamenta
    todos os produtos biologicos veterinarios nos EUA, incluindo
    vacinas. A via para uma nova tecnologia (mRNA) requer:

    Returns:
        Dicionario com passos regulatorios
    """
    return {
        "agency": "USDA-APHIS Center for Veterinary Biologics (CVB)",
        "regulation": "9 CFR Parts 101-118",
        "process": [
            {
                "step": 1,
                "name": "Pre-Submission Meeting",
                "description": (
                    "Reuniao com CVB para discutir nova tecnologia (mRNA-LNP). "
                    "Definir requisitos de seguranca, pureza e potencia."
                ),
                "timeline_months": "1-2",
            },
            {
                "step": 2,
                "name": "Outline of Production",
                "description": (
                    "Documento detalhando processo de producao, controle de "
                    "qualidade e testes de liberacao."
                ),
                "timeline_months": "3-6",
            },
            {
                "step": 3,
                "name": "Experimental Studies (9 CFR 103.3)",
                "description": (
                    "Estudos de seguranca e eficacia em animais alvo (caes). "
                    "Minimo: estudo de seguranca (back-passage, overdose) + "
                    "estudo de eficacia (challenge ou soroconversao)."
                ),
                "timeline_months": "12-24",
            },
            {
                "step": 4,
                "name": "Product License Application",
                "description": (
                    "Submissao de dossiê completo com dados de producao, "
                    "QC, seguranca e eficacia."
                ),
                "timeline_months": "6-12",
            },
            {
                "step": 5,
                "name": "Establishment License",
                "description": (
                    "Licenciamento da planta de producao (GMP veterinario). "
                    "Inspecao in loco pelo CVB."
                ),
                "timeline_months": "6-12",
            },
        ],
        "total_estimated_timeline": "2-5 anos",
        "key_considerations": [
            "Nao existe via acelerada especifica para mRNA veterinario",
            "CVB e receptivo a novas tecnologias (aprovaram vacina de DNA para WNV em equinos)",
            "Dados de seguranca sao mais criticos que eficacia para primeira aprovacao",
            "Possibilidade de conditional license para uso emergencial",
        ],
    }


def _regulatory_brazil() -> dict:
    """Via regulatoria brasileira (MAPA) para biologicos veterinarios.

    O Ministerio da Agricultura regulamenta vacinas veterinarias no
    Brasil atraves da Coordenacao de Fiscalizacao de Produtos
    Veterinarios (CFPV/MAPA).

    Returns:
        Dicionario com requisitos regulatorios brasileiros
    """
    return {
        "agency": "MAPA - Secretaria de Defesa Agropecuaria (SDA)",
        "regulation": "Decreto 5.053/2004 + IN MAPA 44/2007",
        "process": [
            {
                "step": 1,
                "name": "Consulta Previa ao MAPA",
                "description": (
                    "Consulta formal ao Departamento de Saude Animal sobre "
                    "viabilidade regulatoria de vacina de mRNA para caes."
                ),
                "timeline_months": "2-4",
            },
            {
                "step": 2,
                "name": "Desenvolvimento e Testes Pre-Clinicos",
                "description": (
                    "Estudos de seguranca e imunogenicidade em modelo animal. "
                    "Laboratorio credenciado pelo MAPA."
                ),
                "timeline_months": "12-18",
            },
            {
                "step": 3,
                "name": "Submissao de Dossie Tecnico",
                "description": (
                    "Documentacao completa: composicao, processo de producao, "
                    "controle de qualidade, estudos de seguranca e eficacia."
                ),
                "timeline_months": "6-12",
            },
            {
                "step": 4,
                "name": "Avaliacao pela CFPV + Laboratorio Nacional (LANAGRO)",
                "description": (
                    "Analise do dossie e testes laboratoriais independentes. "
                    "LANAGRO verifica identidade, esterilidade, potencia."
                ),
                "timeline_months": "6-18",
            },
            {
                "step": 5,
                "name": "Registro do Produto",
                "description": (
                    "Emissao do registro definitivo apos aprovacao. "
                    "Validade de 5 anos, renovavel."
                ),
                "timeline_months": "3-6",
            },
        ],
        "total_estimated_timeline": "3-5 anos",
        "key_considerations": [
            "Precedente: Leish-Tec (Hertape-Calier) foi registrada em 2014",
            "MAPA nao possui protocolo especifico para mRNA veterinario",
            "Seria necessario criar novo marco regulatorio ou adaptar existente",
            "Parceria com Biomanguinhos facilitaria dialogo com MAPA e ANVISA",
            "Possibilidade de registro experimental para uso em campanhas de saude publica",
        ],
    }


# ---------------------------------------------------------------------------
# Funcao principal
# ---------------------------------------------------------------------------


def run_strategic_analysis() -> StrategicAnalysis:
    """Executa a analise estrategica completa.

    Etapas:
        1. Analisa tamanho de mercado
        2. Constroi argumento de soberania vacinal
        3. Avalia fontes de financiamento e pontua
        4. Mapeia vias regulatorias (Brasil + EUA)
        5. Gera recomendacao estrategica

    Returns:
        StrategicAnalysis com resultados completos
    """
    market = _analyze_market()
    sovereignty = _sovereignty_argument()
    reg_usda = _regulatory_usda()
    reg_brazil = _regulatory_brazil()

    # Formatar fontes de financiamento como dicionarios
    funding_dicts = []
    for fs in FUNDING_SOURCES:
        funding_dicts.append({
            "name": fs.name,
            "institution": fs.institution,
            "description": fs.description,
            "typical_amount_brl_min": fs.typical_amount_brl[0],
            "typical_amount_brl_max": fs.typical_amount_brl[1],
            "probability_pct": fs.probability_pct,
            "timeline_months_min": fs.timeline_months[0],
            "timeline_months_max": fs.timeline_months[1],
            "requirements": fs.requirements,
            "strategic_fit": fs.strategic_fit,
            "score": fs.score,
        })

    # Ranking por score
    ranked = sorted(FUNDING_SOURCES, key=lambda f: f.score, reverse=True)
    ranked_names = [f.name for f in ranked]

    # Recomendacao estrategica
    recommendation = (
        f"Estrategia de financiamento em tres fases:\n"
        f"\n"
        f"FASE 1 (0-12 meses) -- Prova de conceito: ~R$ 500K-2M\n"
        f"  Fonte primaria: {ranked_names[2]} (CNPq/MCTI)\n"
        f"  Fonte secundaria: {ranked_names[3]} (FAP estadual)\n"
        f"  Objetivo: dados pre-clinicos em camundongos BALB/c\n"
        f"\n"
        f"FASE 2 (12-30 meses) -- Pre-clinico avancado: ~R$ 5-10M\n"
        f"  Fonte primaria: {ranked_names[1]} (FINEP subvencao)\n"
        f"  Fonte secundaria: {ranked_names[0]} (Fiocruz PDT)\n"
        f"  Objetivo: eficacia em modelo canino, GLP safety\n"
        f"\n"
        f"FASE 3 (30-60 meses) -- Registro e escala: ~R$ 20-50M\n"
        f"  Fonte primaria: {ranked_names[0]} (Fiocruz/Biomanguinhos)\n"
        f"  Fonte secundaria: BNDES Saude Animal\n"
        f"  Objetivo: registro MAPA, transferencia de tecnologia\n"
        f"\n"
        f"Mercado total enderecavel: ~{market.addressable_market_dogs:,} caes/ano "
        f"(R$ {market.market_value_brl_year:,.0f}/ano). "
        f"O argumento One Health (vacinar caes para proteger humanos) "
        f"e o diferencial estrategico para justificar investimento publico."
    )

    return StrategicAnalysis(
        market=market,
        sovereignty=sovereignty,
        funding_sources=funding_dicts,
        funding_ranked=ranked_names,
        regulatory_usda=reg_usda,
        regulatory_brazil=reg_brazil,
        recommendation=recommendation,
    )


def to_dict(analysis: StrategicAnalysis) -> dict:
    """Serializa a analise estrategica para JSON.

    Args:
        analysis: resultado da analise

    Returns:
        Dicionario serializavel para JSON
    """
    m = analysis.market

    return {
        "market": {
            "total_dogs_brazil": m.total_dogs_brazil,
            "dogs_at_risk": m.dogs_at_risk,
            "dogs_at_risk_fraction_pct": round(m.dogs_at_risk_fraction * 100, 1),
            "seroprevalence_mean_pct": round(m.seroprevalence_mean * 100, 1),
            "estimated_infected_dogs": m.estimated_infected_dogs,
            "addressable_market_dogs": m.addressable_market_dogs,
            "addressable_market_doses_year": m.addressable_market_doses,
            "market_value_brl_year": m.market_value_brl_year,
            "market_value_usd_year": m.market_value_usd_year,
            "human_cases_per_year": m.human_cases_per_year,
            "human_fatality_rate_pct": round(m.human_fatality_rate * 100, 1),
            "annual_control_cost_brl": m.annual_control_cost_brl,
            "competitor_leish_tec_price_brl": m.competitor_price_brl,
        },
        "sovereignty_argument": {
            "title": analysis.sovereignty.title,
            "context": analysis.sovereignty.context,
            "covid19_precedent": analysis.sovereignty.precedent,
            "brazil_mrna_capacity": analysis.sovereignty.brazil_mrna_capacity,
            "one_health": analysis.sovereignty.one_health_argument,
            "strategic_value": analysis.sovereignty.strategic_value,
        },
        "funding_sources": analysis.funding_sources,
        "funding_ranked_by_score": analysis.funding_ranked,
        "regulatory_usda": analysis.regulatory_usda,
        "regulatory_brazil": analysis.regulatory_brazil,
        "recommendation": analysis.recommendation,
    }
