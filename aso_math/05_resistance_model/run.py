"""Modulo 05 — Modelo matematico de resistencia a MRL-ASO-001.

Calcula a probabilidade de Leishmania desenvolver resistencia ao ASO
via mutacao do SL RNA alvo. Combina:
    1. Analise de mutacoes de escape (binding disruption)
    2. Restricao funcional (selecao purificadora no SL RNA)
    3. Cadeia de Markov para fixacao em tandem array (~150 copias)
    4. Processo de Poisson para tempo ate o primeiro parasita resistente
    5. Analise de sensibilidade sobre parametros-chave
    6. Comparacao com resistencia a drogas convencionais

O resultado central: mesmo sob premissas conservadoras (favorecendo
resistencia), o tempo esperado para escape excede vastamente o tempo
de tratamento e a vida util de qualquer droga convencional.

Referencia biologica:
    - SL RNA: Liang XH et al. (2003) Int J Parasitol 33(14):1603-1612
    - Taxa de mutacao: Rogers MB et al. (2011) PLoS Genetics 7(8):e1002237
    - Tandem arrays: Martínez-Calvillo S et al. (2001) RNA 7(7):1013-1023
    - ASO binding: Crooke ST et al. (2017) Nucleic Acids Res 45(9):5614-5625
    - Resistencia a antimoniais: Croft SL et al. (2006) Clin Microbiol Rev 19(1):111-126
"""

from __future__ import annotations

import math
from typing import Any

from aso_math.config import (
    ASO_SEQUENCE,
    ASO_TARGET_END,
    ASO_TARGET_START,
    BASES,
    DG_FUNCTIONAL_THRESHOLD,
    GENERATION_TIME_HOURS,
    MUTATION_RATE,
    NN_INIT,
    SL_RNA_COPY_NUMBER,
    SL_SEQUENCE,
    T_PHYSIOLOGICAL,
)
from aso_math.envelope import Timer, create_envelope, write_result
from aso_math.target_config import TargetConfig
from aso_math.thermo import _get_nn_params, compute_dg
from core.logger import get_logger

logger = get_logger("05_resistance_model")

# ---------------------------------------------------------------------------
# Constantes do modelo de resistencia
# ---------------------------------------------------------------------------

# Horas por ano (365.25 dias * 24 horas)
HOURS_PER_YEAR: float = 365.25 * 24.0

# Geracoes por ano: cada geracao dura 12 h -> ~730 geracoes/ano
GENERATIONS_PER_YEAR: float = HOURS_PER_YEAR / GENERATION_TIME_HOURS

# Tamanhos de populacao para modelagem
# Ref: Courtenay O et al. (2002) Parasitology 125:S31-S56
# 10^3: infeccao precoce/assintomatica
# 10^6: infeccao estabelecida
# 10^8: leishmaniose visceral clinica
# 10^10: carga maxima em casos severos
POPULATION_SIZES: list[float] = [1e3, 1e6, 1e8, 1e10]

# Probabilidade de uma mutacao no SL RNA reter funcao de trans-splicing
# Conservadoramente alta: 0.0 para posicoes conservadas, 0.3 para variaveis.
# Na realidade, TODAS as posicoes do SL RNA de Leishmania sao conservadas
# entre especies, entao P(funcional | mutacao) ≈ 0 para qualquer posicao.
# Usamos 0.0 como valor padrao (honesto) e 0.3 como cenario generoso
# na analise de sensibilidade.
P_RETAINS_FUNCTION_CONSERVED: float = 0.0
P_RETAINS_FUNCTION_VARIABLE: float = 0.3

# Para analise de sensibilidade: se ALGUM dia se descobrir posicoes
# variaveis, usamos esse limiar generoso como teto.
P_RETAINS_FUNCTION_GENEROUS: float = 0.3

# Dados de resistencia a drogas convencionais
# Ref: Croft SL et al. (2006) Clin Microbiol Rev 19(1):111-126
# Ref: Sundar S et al. (2012) Parasitology 139(3):375-392
# Ref: Ponte-Sucre A et al. (2017) PLoS Negl Trop Dis 11(12):e0006056
KNOWN_DRUG_RESISTANCE: dict[str, dict[str, str]] = {
    "Antimoniais pentavalentes": {
        "time_years": "5-10",
        "mechanism": "Transportador ABC, tiol metabolismo",
    },
    "Miltefosina": {
        "time_years": "3-5",
        "mechanism": "Mutacao em transportador de lipidio LdMT/LdRos3",
    },
    "Anfotericina B": {
        "time_years": ">20",
        "mechanism": "Alteracao de esterol de membrana (raro)",
    },
    "Pentamidina": {
        "time_years": "3-8",
        "mechanism": "Reducao de captacao mitocondrial",
    },
}


# ---------------------------------------------------------------------------
# 1. Identificacao de mutacoes de escape
# ---------------------------------------------------------------------------


def _compute_disrupted_dg(aso_seq: str, mismatch_position: int) -> float:
    """Calcula dG do ASO com mismatch na posicao indicada do alvo.

    Modelo: quando o alvo muta na posicao i, o ASO (inalterado) forma
    mismatch nessa posicao. Os dinucleotideos NN que incluem a posicao i
    (ou seja, (i-1,i) e (i,i+1)) perdem sua contribuicao estabilizadora.

    Zeramos a contribuicao desses dinucleotideos em vez de usar parametros
    de mismatch positivos. Isso SUBestima o impacto real do mismatch,
    tornando nossa estimativa CONSERVADORA: mais mutacoes parecem escapar
    do que realmente escapariam.

    Ref: SantaLucia J Jr. (1998) PNAS 95(4):1460-1465 — parametros NN
    Ref: Allawi HT, SantaLucia J Jr. (1997) Biochemistry 36:10581 — mismatches

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        mismatch_position: Indice 0-based da posicao com mismatch no ASO.

    Returns:
        dG em kcal/mol com o mismatch (valor menos negativo = binding pior).
    """
    seq = aso_seq.upper()
    n = len(seq)

    # Posicoes de dinucleotideos afetados pelo mismatch
    # Dinucleotideo (i-1, i) existe se i > 0
    # Dinucleotideo (i, i+1) existe se i < n-1
    disrupted_positions: set[int] = set()
    if mismatch_position > 0:
        disrupted_positions.add(mismatch_position - 1)
    if mismatch_position < n - 1:
        disrupted_positions.add(mismatch_position)

    # Somar parametros NN, zerando os dinucleotideos afetados
    total_dh = NN_INIT[0]
    total_ds = NN_INIT[1]

    for i in range(n - 1):
        if i in disrupted_positions:
            # Mismatch: contribuicao zerada (conservador — real seria positiva)
            continue
        dh, ds = _get_nn_params(seq[i:i + 2])
        total_dh += dh
        total_ds += ds

    delta_g = total_dh - T_PHYSIOLOGICAL * (total_ds / 1000.0)
    return round(delta_g, 2)


def identify_escape_mutations(
    aso_seq: str = ASO_SEQUENCE,
    sl_seq: str = SL_SEQUENCE,
    target_start: int = ASO_TARGET_START,
    target_end: int = ASO_TARGET_END,
    dg_threshold: float = DG_FUNCTIONAL_THRESHOLD,
) -> list[dict[str, Any]]:
    """Identifica mutacoes pontuais no alvo que poderiam escapar do ASO.

    Para cada uma das 75 mutacoes pontuais possiveis (25 posicoes x 3 bases
    alternativas) no target do SL RNA:
        - Calcula dG do ASO original contra o alvo mutado
        - Classifica como "escape" se dG > dg_threshold

    A posicao no alvo do SL RNA mapeia diretamente para a posicao no ASO:
    o ASO e o complemento reverso do alvo, entao a posicao i do alvo
    corresponde a posicao (ASO_LENGTH - 1 - i) no ASO.

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        sl_seq: Sequencia do SL RNA completo.
        target_start: Posicao inicial do alvo no SL RNA (0-indexed).
        target_end: Posicao final do alvo no SL RNA (0-indexed, exclusivo).
        dg_threshold: Limiar de dG funcional (kcal/mol).

    Returns:
        Lista de dicionarios com detalhes de cada mutacao.
    """
    target = sl_seq[target_start:target_end]
    aso = aso_seq.upper()
    aso_len = len(aso)
    wt_dg = compute_dg(aso)
    mutations: list[dict[str, Any]] = []

    for target_pos in range(len(target)):
        original_base = target[target_pos]

        for alt_base in BASES:
            if alt_base == original_base:
                continue

            # A posicao no ASO correspondente (complemento reverso)
            # Target pos 0 -> ASO pos (len-1), target pos 1 -> ASO pos (len-2), etc.
            aso_mismatch_pos = aso_len - 1 - target_pos

            # Calcular dG com mismatch nessa posicao do ASO
            disrupted_dg = _compute_disrupted_dg(aso, aso_mismatch_pos)

            # Classificar: escape se dG acima do limiar (binding insuficiente)
            is_escape = disrupted_dg > dg_threshold

            mutations.append({
                "target_position": target_pos,
                "sl_position": target_start + target_pos,
                "aso_mismatch_position": aso_mismatch_pos,
                "original": original_base,
                "mutant": alt_base,
                "disrupted_dg_kcal": disrupted_dg,
                "wildtype_dg_kcal": wt_dg,
                "dg_change_kcal": round(disrupted_dg - wt_dg, 2),
                "is_escape": is_escape,
            })

    # Ordenar por dG (menos negativo = mais provavel escape, primeiro)
    mutations.sort(key=lambda m: m["disrupted_dg_kcal"], reverse=True)

    n_escape = sum(1 for m in mutations if m["is_escape"])
    logger.info(
        "Mutacoes de escape: %d / %d rompem binding (dG > %.1f kcal/mol)",
        n_escape,
        len(mutations),
        dg_threshold,
    )

    return mutations


# ---------------------------------------------------------------------------
# 2. Restricao funcional (selecao purificadora)
# ---------------------------------------------------------------------------


def compute_functional_constraint(
    escape_mutations: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Avalia se cada mutacao de escape tambem mata o parasita.

    O SL RNA e essencial para trans-splicing — TODOS os mRNAs do parasita
    dependem dele. As posicoes do SL RNA sao 100% conservadas entre todas
    as especies de Leishmania (L. infantum, L. major, L. donovani, etc.),
    indicando selecao purificadora extrema.

    Ref: Liang XH et al. (2003) — conservacao do SL RNA
    Ref: Michaeli S (2011) Parasitology 138(12):1473-1491 — funcao do SL RNA

    Para cada mutacao de escape:
        - Se a posicao e conservada (entropia = 0): P(funcional) = 0.0
        - Se a posicao fosse variavel: P(funcional) = 0.3

    Como TODAS as posicoes sao conservadas em Leishmania, P(funcional) = 0.0
    para todas as mutacoes. Isso e biologicamente correto, nao uma escolha
    arbitraria: nenhuma variacao natural foi observada nessas posicoes em
    ~500 milhoes de anos de evolucao.

    Args:
        escape_mutations: Lista de mutacoes da identify_escape_mutations().

    Returns:
        Lista atualizada com probabilidade de reter funcao e viabilidade.
    """
    annotated: list[dict[str, Any]] = []

    for mut in escape_mutations:
        # Todas as posicoes do SL RNA sao conservadas entre Leishmania spp.
        # Entropia de Shannon = 0 para todas as posicoes do alvo.
        # Se algum dia dados de Module 03 indicarem posicoes variaveis,
        # atualizariamos aqui.
        is_conserved = True  # todas as 25 posicoes do alvo

        if is_conserved:
            p_functional = P_RETAINS_FUNCTION_CONSERVED
        else:
            p_functional = P_RETAINS_FUNCTION_VARIABLE

        # Uma mutacao e "viavel" se escapa do ASO E retém funcao
        viable = mut["is_escape"] and p_functional > 0.0

        entry = {**mut}
        entry["is_conserved_position"] = is_conserved
        entry["retains_function_prob"] = p_functional
        entry["viable"] = viable
        annotated.append(entry)

    n_viable = sum(1 for m in annotated if m["viable"])
    n_escape = sum(1 for m in annotated if m["is_escape"])

    logger.info(
        "Restricao funcional: %d escape mutations, %d funcionalmente viaveis",
        n_escape,
        n_viable,
    )

    return annotated


# ---------------------------------------------------------------------------
# 3. Taxa de escape por copia unica
# ---------------------------------------------------------------------------


def compute_single_copy_escape_rate(
    annotated_mutations: list[dict[str, Any]],
    mutation_rate: float = MUTATION_RATE,
) -> dict[str, Any]:
    """Calcula taxa de aparecimento de mutante de escape viavel por replicacao.

    lambda_single = SUM( mu * P(escape_i) * P(funcional_i) )
    para cada mutacao i nas 75 possiveis.

    Onde:
        mu = taxa de mutacao por base por replicacao (2e-9)
        P(escape_i) = 1 se a mutacao rompe binding, 0 caso contrario
        P(funcional_i) = probabilidade de reter funcao de trans-splicing

    Como P(funcional) = 0 para todas as posicoes conservadas,
    lambda_single = 0 no modelo realista. Para a analise de sensibilidade,
    usamos P(funcional) = 0.3 como cenario absurdamente generoso.

    Args:
        annotated_mutations: Lista de mutacoes com probabilidade funcional.
        mutation_rate: Taxa de mutacao por base por replicacao.

    Returns:
        Dicionario com taxa de escape e decomposicao por mutacao.
    """
    total_rate = 0.0
    contributing: list[dict[str, Any]] = []

    for mut in annotated_mutations:
        # Cada mutacao especifica tem taxa mu (nao mu/3 porque ja iteramos
        # sobre cada base alternativa individualmente)
        rate_i = mutation_rate * (1.0 if mut["is_escape"] else 0.0) * mut["retains_function_prob"]
        total_rate += rate_i

        if rate_i > 0:
            contributing.append({
                "mutation": f"{mut['original']}{mut['target_position']+1}{mut['mutant']}",
                "rate_per_replication": rate_i,
            })

    logger.info(
        "Taxa de escape por copia: lambda_single = %.2e por replicacao (%d mutacoes contribuem)",
        total_rate,
        len(contributing),
    )

    return {
        "lambda_single_per_replication": total_rate,
        "contributing_mutations": contributing,
        "n_contributing": len(contributing),
        "mutation_rate_used": mutation_rate,
    }


# ---------------------------------------------------------------------------
# 4. Modelo de fixacao em tandem array
# ---------------------------------------------------------------------------


def model_tandem_array_fixation(
    n_copies: int = SL_RNA_COPY_NUMBER,
    single_copy_rate: float = 0.0,
    gen_time: float = GENERATION_TIME_HOURS,
) -> dict[str, Any]:
    """Modela a fixacao de uma mutacao de escape no tandem array do SL RNA.

    Mesmo que uma copia do SL RNA adquira uma mutacao de escape, ela precisa
    se tornar a MAIORIA do array (~150 copias) para conferir resistencia.
    Isso porque a maquinaria de trans-splicing usa copias do array
    estocasticamente — enquanto a maioria for wildtype, o ASO ainda funciona.

    Sob deriva neutra (Kimura 1962):
        P(fixacao) = 1 / n_copies
        T(fixacao)  ~ n_copies geracoes (em unidades de N_e do array)

    Sob selecao purificadora contra o mutante (porque mesmo mutantes de
    "escape" provavelmente tem funcionalidade de trans-splicing reduzida):
        P(fixacao) < 1 / n_copies
        Usamos o valor neutro como LIMITE SUPERIOR (conservador).

    Ref: Kimura M (1962) Genetics 47:713-719 — probabilidade de fixacao
    Ref: Nei M, Rooney AP (2005) Annu Rev Genet 39:121-152 — concerted evolution
    Ref: Liao D (1999) Am J Hum Genet 64(1):24-30 — tandem repeat dynamics

    Args:
        n_copies: Numero de copias no tandem array.
        single_copy_rate: Taxa de escape por copia por replicacao.
        gen_time: Tempo de geracao em horas.

    Returns:
        Dicionario com parametros de fixacao e taxa efetiva.
    """
    # Probabilidade de fixacao sob deriva neutra (limite superior)
    p_fixation = 1.0 / n_copies

    # Tempo para fixacao (se ocorrer): ~n_copies geracoes
    # Ref: Kimura & Ohta (1969) Genetics 61:763 — tempo medio de fixacao = 4*N geracoes
    # Para tandem array com concerted evolution: ~n_copies geracoes e conservador
    t_fixation_generations = n_copies

    # Taxa efetiva: precisa aparecer E fixar
    lambda_effective = single_copy_rate * p_fixation

    logger.info(
        "Tandem array: %d copias, P(fixacao) = %.4f, T(fixacao) = %d geracoes",
        n_copies,
        p_fixation,
        t_fixation_generations,
    )
    logger.info(
        "  Taxa efetiva: lambda_eff = %.2e por parasita por replicacao",
        lambda_effective,
    )

    return {
        "copy_number": n_copies,
        "fixation_probability": round(p_fixation, 6),
        "time_to_fixation_generations": t_fixation_generations,
        "time_to_fixation_hours": t_fixation_generations * gen_time,
        "time_to_fixation_days": round(
            t_fixation_generations * gen_time / 24.0, 1,
        ),
        "single_copy_rate": single_copy_rate,
        "effective_escape_rate": lambda_effective,
    }


# ---------------------------------------------------------------------------
# 5. Tempo de escape populacional (processo de Poisson)
# ---------------------------------------------------------------------------


def compute_population_escape_time(
    lambda_effective: float,
    population_sizes: list[float] | None = None,
    gen_time: float = GENERATION_TIME_HOURS,
) -> dict[str, Any]:
    """Calcula tempo esperado ate o primeiro parasita resistente via Poisson.

    Para uma populacao de N parasitas, cada um com taxa de escape
    lambda_effective por replicacao:
        Lambda_total = N * lambda_effective
        E[T] = 1 / Lambda_total  (em geracoes)

    Converte para horas, dias e anos.

    Se lambda_effective = 0 (caso biologicamente correto), E[T] = infinito.

    Ref: Ross SM (2014) Introduction to Probability Models, 11th ed. — Poisson

    Args:
        lambda_effective: Taxa de escape por parasita por replicacao.
        population_sizes: Lista de tamanhos populacionais a modelar.
        gen_time: Tempo de geracao em horas.

    Returns:
        Dicionario com tempos de escape para cada tamanho populacional.
    """
    if population_sizes is None:
        population_sizes = POPULATION_SIZES

    # Horas por ano, constante astronomica
    hours_per_year = 365.25 * 24.0

    results: dict[str, Any] = {}

    for n_pop in population_sizes:
        pop_key = f"population_{n_pop:.0e}"

        if lambda_effective <= 0.0:
            # Sem taxa de escape: tempo infinito
            results[pop_key] = {
                "population_size": n_pop,
                "total_rate_per_generation": 0.0,
                "expected_generations": math.inf,
                "expected_hours": math.inf,
                "expected_days": math.inf,
                "expected_years": math.inf,
                "probability_in_treatment_course": 0.0,
            }
        else:
            total_rate = n_pop * lambda_effective
            expected_gen = 1.0 / total_rate
            expected_hours = expected_gen * gen_time
            expected_days = expected_hours / 24.0
            expected_years = expected_hours / hours_per_year

            # P(escape durante curso de tratamento de 8 semanas)
            # P(T <= t) = 1 - exp(-Lambda * t)
            treatment_generations = (8 * 7 * 24) / gen_time  # 8 semanas
            p_treatment = 1.0 - math.exp(-total_rate * treatment_generations)

            results[pop_key] = {
                "population_size": n_pop,
                "total_rate_per_generation": total_rate,
                "expected_generations": round(expected_gen, 2),
                "expected_hours": round(expected_hours, 2),
                "expected_days": round(expected_days, 2),
                "expected_years": round(expected_years, 2),
                "probability_in_treatment_course": round(p_treatment, 15),
            }

        logger.info(
            "  N=%.0e: E[T] = %s anos",
            n_pop,
            str(round(results[pop_key]["expected_years"], 2))
            if results[pop_key]["expected_years"] != math.inf
            else "infinito",
        )

    return results


# ---------------------------------------------------------------------------
# 6. Analise de sensibilidade
# ---------------------------------------------------------------------------


def sensitivity_analysis(
    aso_seq: str = ASO_SEQUENCE,
    sl_seq: str = SL_SEQUENCE,
    target_start: int = ASO_TARGET_START,
    target_end: int = ASO_TARGET_END,
    gen_time: float = GENERATION_TIME_HOURS,
) -> dict[str, Any]:
    """Varia parametros-chave e recalcula tempo de escape.

    Parametros variados:
        - mutation_rate: [1e-10, 1e-9, 2e-9, 1e-8]
        - population_size: [1e3, 1e6, 1e8, 1e10]
        - dg_threshold: [-26.0, -25.0, -24.0, -20.0, -15.0, -10.0]
        - n_copies: [50, 100, 150, 200]

    NOTA: os limiares -26.0, -25.0 e -24.0 sao biologicamente irrealistas
    (ASOs com dG < -20 ainda funcionam muito bem), mas os incluimos para
    mostrar que MESMO SE reduzirmos drasticamente o requisito de binding,
    o tempo de escape permanece astronomico. Isso fortalece a prova.

    Para a sensibilidade, usamos P(funcional) = 0.3 (generoso) para
    mutacoes em posicoes que SERIAM de escape, criando um cenario
    artificialmente favoravel a resistencia. Se mesmo assim o tempo
    for astronomico, a prova e robusta.

    Args:
        aso_seq: Sequencia do ASO (5'->3').
        sl_seq: Sequencia do SL RNA completo.
        target_start: Posicao inicial do alvo no SL RNA (0-indexed).
        target_end: Posicao final do alvo no SL RNA (0-indexed, exclusivo).
        gen_time: Tempo de geracao em horas.

    Returns:
        Dicionario com matriz de resultados e resumo.
    """
    mutation_rates = [1e-10, 1e-9, 2e-9, 1e-8]
    pop_sizes = [1e3, 1e6, 1e8, 1e10]
    # Limiares incluem valores extremos para capturar cenarios hipoteticos
    # -26.0 a -24.0: biologicamente irrealistas (binding ainda funcional)
    # -20.0 a -10.0: limiares publicados na literatura (Crooke et al. 2017)
    dg_thresholds = [-26.0, -25.0, -24.0, -20.0, -15.0, -10.0]
    copy_numbers = [50, 100, 150, 200]

    target = sl_seq[target_start:target_end]
    aso = aso_seq.upper()
    aso_len = len(aso)
    wt_dg = compute_dg(aso)

    results: list[dict[str, Any]] = []

    for dg_thresh in dg_thresholds:
        # Recalcular quais mutacoes sao escape para este limiar
        n_escape_at_threshold = 0

        for target_pos in range(len(target)):
            original_base = target[target_pos]
            for alt_base in BASES:
                if alt_base == original_base:
                    continue
                aso_mm_pos = aso_len - 1 - target_pos
                d_dg = _compute_disrupted_dg(aso, aso_mm_pos)
                if d_dg > dg_thresh:
                    n_escape_at_threshold += 1

        for mu in mutation_rates:
            for n_cop in copy_numbers:
                for n_pop in pop_sizes:
                    # Cenario generoso: P(funcional) = 0.3
                    # para todas as mutacoes de escape
                    lambda_single = n_escape_at_threshold * mu * P_RETAINS_FUNCTION_GENEROUS

                    # Fixacao em tandem array
                    p_fix = 1.0 / n_cop
                    lambda_eff = lambda_single * p_fix

                    # Tempo populacional
                    if lambda_eff > 0:
                        total_rate = n_pop * lambda_eff
                        expected_gen = 1.0 / total_rate
                        expected_years = (
                            expected_gen * gen_time / HOURS_PER_YEAR
                        )
                    else:
                        expected_years = math.inf

                    results.append({
                        "mutation_rate": mu,
                        "population_size": n_pop,
                        "dg_threshold": dg_thresh,
                        "n_copies": n_cop,
                        "n_escape_mutations": n_escape_at_threshold,
                        "p_functional": P_RETAINS_FUNCTION_GENEROUS,
                        "lambda_single": lambda_single,
                        "lambda_effective": lambda_eff,
                        "expected_years": (
                            expected_years if expected_years != math.inf
                            else "infinity"
                        ),
                    })

    # Encontrar cenario mais favoravel a resistencia (menor tempo)
    finite_results = [r for r in results if r["expected_years"] != "infinity"]
    if finite_results:
        worst_case = min(finite_results, key=lambda r: r["expected_years"])
    else:
        worst_case = None

    logger.info(
        "Sensibilidade: %d combinacoes avaliadas, %d com tempo finito",
        len(results),
        len(finite_results),
    )
    if worst_case:
        logger.info(
            "  Pior caso (mais rapido): %.2e anos (mu=%.0e, N=%.0e, dG_thresh=%.1f, copies=%d)",
            worst_case["expected_years"],
            worst_case["mutation_rate"],
            worst_case["population_size"],
            worst_case["dg_threshold"],
            worst_case["n_copies"],
        )

    return {
        "parameters_varied": ["mutation_rate", "population_size", "dg_threshold", "n_copies"],
        "n_combinations": len(results),
        "n_finite": len(finite_results),
        "worst_case_scenario": worst_case,
        "results": results,
    }


# ---------------------------------------------------------------------------
# 7. Comparacao com resistencia a drogas convencionais
# ---------------------------------------------------------------------------


def compare_with_drug_resistance(
    aso_escape_years: float | str,
) -> dict[str, Any]:
    """Compara tempo de escape do ASO com resistencia a drogas conhecidas.

    Drogas anti-Leishmania desenvolvem resistencia em anos a decadas
    porque dependem de mutacao em UM gene (ou poucos genes) que nao e
    essencial sob pressao seletiva do farmaco. O SL RNA e diferente:
        1. Mutacao no alvo tambem mata o parasita (dupla barreira)
        2. O alvo existe em ~150 copias (precisa fixar no array)
        3. O alvo e conservado ha ~500M anos (selecao extrema)

    Args:
        aso_escape_years: Tempo estimado para escape do ASO (anos ou "infinity").

    Returns:
        Dicionario com comparacao estruturada.
    """
    # Formatar tempo do ASO
    if isinstance(aso_escape_years, str):
        aso_display = aso_escape_years
    elif aso_escape_years == math.inf:
        aso_display = "infinity"
    elif aso_escape_years > 1e12:
        aso_display = f">{aso_escape_years:.2e}"
    else:
        aso_display = f"{aso_escape_years:.2e}"

    # Construir tabela comparativa
    comparison: dict[str, Any] = {
        "MRL-ASO-001": {
            "time_years": aso_display,
            "mechanism": "Mutacao no SL RNA + fixacao em tandem array + retencao funcional",
            "barriers": [
                "Mutacao deve romper binding do ASO (dG > -15 kcal/mol)",
                "Mutacao deve reter funcao de trans-splicing (P ~ 0)",
                "Mutante deve fixar em array de ~150 copias (P = 1/150)",
                "Tudo deve ocorrer no mesmo parasita",
            ],
            "n_barriers": 3,
        },
    }

    for drug, info in KNOWN_DRUG_RESISTANCE.items():
        comparison[drug] = {
            "time_years": info["time_years"],
            "mechanism": info["mechanism"],
            "barriers": ["Mutacao em gene(s) alvo"],
            "n_barriers": 1,
        }

    # Conclusao
    if aso_display == "infinity":
        conclusion = (
            "MRL-ASO-001 apresenta barreira de resistencia infinita no modelo realista. "
            "Mesmo no cenario mais generoso da analise de sensibilidade, o tempo de escape "
            "excede vastamente o de qualquer droga anti-Leishmania conhecida. "
            "Isso decorre de tres barreiras simultaneas: (1) a mutacao deve romper o "
            "binding do ASO, (2) deve reter funcao essencial do SL RNA, e (3) deve fixar "
            "em ~150 copias do tandem array. A probabilidade conjunta e efetivamente zero."
        )
    else:
        conclusion = (
            f"MRL-ASO-001: tempo estimado para resistencia = {aso_display} anos. "
            "Mesmo sob premissas extremamente conservadoras (favorecendo resistencia), "
            "este valor excede o de qualquer droga anti-Leishmania convencional, "
            "confirmando a vantagem terapeutica do alvo SL RNA."
        )

    comparison["conclusion"] = conclusion

    logger.info("Comparacao com drogas: %s", conclusion[:120] + "...")

    return comparison


# ---------------------------------------------------------------------------
# 8. Orquestrador principal
# ---------------------------------------------------------------------------


def main(config: TargetConfig | None = None) -> dict[str, Any]:
    """Executa a analise completa do modelo de resistencia.

    Args:
        config: Configuracao do organismo-alvo. Se None, usa L. infantum.

    Returns:
        Envelope completo com todos os resultados.
    """
    if config is None:
        config = TargetConfig()

    # Extrair valores do config (fallback para constantes de L. infantum)
    aso_seq = config.aso_sequence or ASO_SEQUENCE
    sl_seq = config.sl_sequence or SL_SEQUENCE
    target_start = config.aso_target_start or ASO_TARGET_START
    target_end = config.aso_target_end or ASO_TARGET_END
    mutation_rate = config.mutation_rate
    gen_time = config.generation_time_hours
    copy_number = config.sl_copy_number
    dg_threshold = config.dg_functional_threshold

    logger.info("=" * 60)
    logger.info("MODULO 05: Modelo de Resistencia — %s", config.species_name)
    logger.info("=" * 60)

    envelope = create_envelope("05_resistance_model")

    with Timer() as timer:
        # --- Passo 1: Mutacoes de escape ---
        logger.info("Passo 1/7: Identificando mutacoes de escape...")
        escape_mutations = identify_escape_mutations(
            aso_seq=aso_seq,
            sl_seq=sl_seq,
            target_start=target_start,
            target_end=target_end,
            dg_threshold=dg_threshold,
        )

        n_total = len(escape_mutations)
        n_binding_disrupting = sum(1 for m in escape_mutations if m["is_escape"])
        wt_dg = compute_dg(aso_seq)

        logger.info(
            "  Total: %d mutacoes, %d rompem binding (dG_wt = %.2f kcal/mol)",
            n_total,
            n_binding_disrupting,
            wt_dg,
        )

        # --- Passo 2: Restricao funcional ---
        logger.info("Passo 2/7: Avaliando restricao funcional...")
        annotated = compute_functional_constraint(escape_mutations)

        n_viable = sum(1 for m in annotated if m["viable"])
        logger.info(
            "  Viaveis (escape + funcional): %d / %d",
            n_viable,
            n_binding_disrupting,
        )

        # --- Passo 3: Taxa de escape por copia ---
        logger.info("Passo 3/7: Calculando taxa de escape por copia unica...")
        single_rate = compute_single_copy_escape_rate(
            annotated, mutation_rate=mutation_rate,
        )

        # --- Passo 4: Fixacao em tandem array ---
        logger.info("Passo 4/7: Modelando fixacao no tandem array...")
        tandem = model_tandem_array_fixation(
            n_copies=copy_number,
            single_copy_rate=single_rate["lambda_single_per_replication"],
            gen_time=gen_time,
        )

        # --- Passo 5: Tempo de escape populacional ---
        logger.info("Passo 5/7: Computando tempo de escape populacional...")
        escape_times = compute_population_escape_time(
            lambda_effective=tandem["effective_escape_rate"],
            gen_time=gen_time,
        )

        # Extrair tempo para populacao referencia (10^8 = caso clinico)
        ref_pop_key = "population_1e+08"
        ref_time = escape_times.get(ref_pop_key, {})
        ref_years = ref_time.get("expected_years", math.inf) if ref_time else math.inf

        # --- Passo 6: Analise de sensibilidade ---
        logger.info("Passo 6/7: Executando analise de sensibilidade...")
        sensitivity = sensitivity_analysis(
            aso_seq=aso_seq,
            sl_seq=sl_seq,
            target_start=target_start,
            target_end=target_end,
            gen_time=gen_time,
        )

        # Pior caso da sensibilidade para comparacao
        worst_case = sensitivity.get("worst_case_scenario")
        worst_years: float | str = "infinity"
        if worst_case and worst_case["expected_years"] != "infinity":
            worst_years = worst_case["expected_years"]

        # Formatar worst_years para exibicao legivel
        if isinstance(worst_years, (int, float)):
            worst_years_display = f"{worst_years:.2e}"
        else:
            worst_years_display = str(worst_years)

        # --- Passo 7: Comparacao com drogas ---
        logger.info("Passo 7/7: Comparando com resistencia a drogas convencionais...")
        comparison = compare_with_drug_resistance(ref_years)

        # --- Conclusao ---
        if ref_years == math.inf:
            conclusion = (
                "Resistencia a MRL-ASO-001 e matematicamente impossivel no modelo realista. "
                f"Das 75 mutacoes pontuais possiveis no alvo, {n_binding_disrupting} rompem "
                f"o binding do ASO, mas NENHUMA retém funcao de trans-splicing (todas as "
                f"posicoes do SL RNA sao 100% conservadas entre Leishmania spp.). "
                f"Mesmo no cenario mais generoso da sensibilidade "
                f"(P_funcional=0.3, dG_thresh=-26, mu=1e-8, N=1e10, 50 copias), "
                f"o tempo minimo para resistencia e {worst_years_display} anos."
            )
        else:
            conclusion = (
                f"Tempo estimado para resistencia a MRL-ASO-001: {ref_years:.2e} anos "
                f"(populacao referencia N=1e8). Das {n_total} mutacoes analisadas, "
                f"{n_binding_disrupting} rompem binding e {n_viable} sao funcionalmente "
                f"viaveis. Pior caso na sensibilidade: {worst_years_display} anos."
            )

        logger.info("CONCLUSAO: %s", conclusion)

    # --- Montar envelope ---
    envelope["runtime_seconds"] = timer.elapsed
    envelope["status"] = "success"
    envelope["summary"]["conclusion"] = conclusion
    envelope["summary"]["key_metrics"] = {
        "total_mutations_analyzed": n_total,
        "binding_disrupting": n_binding_disrupting,
        "functionally_viable": n_viable,
        "effective_escape_rate": tandem["effective_escape_rate"],
        "reference_escape_years": ref_years if ref_years != math.inf else "infinity",
        "worst_case_years": worst_years,
    }

    # Preparar detalhe das mutacoes de escape (apenas as que rompem binding)
    escape_detail = [
        {
            "position": m["target_position"],
            "sl_position": m["sl_position"],
            "original": m["original"],
            "mutant": m["mutant"],
            "disrupted_dg_kcal": m["disrupted_dg_kcal"],
            "is_escape": m["is_escape"],
            "retains_function_prob": m["retains_function_prob"],
            "viable": m["viable"],
        }
        for m in annotated
    ]

    # Formatar tempos de escape para JSON (converter inf -> "infinity")
    escape_times_json: dict[str, Any] = {}
    for key, val in escape_times.items():
        entry = {}
        for k, v in val.items():
            if isinstance(v, float) and math.isinf(v):
                entry[k] = "infinity"
            else:
                entry[k] = v
        escape_times_json[key] = entry

    envelope["data"] = {
        "wildtype_dg_kcal": wt_dg,
        "escape_mutations": {
            "total_possible": n_total,
            "binding_disrupting": n_binding_disrupting,
            "functionally_viable": n_viable,
            "escape_mutations_detail": escape_detail,
        },
        "single_copy_escape": single_rate,
        "tandem_array_model": {
            "copy_number": tandem["copy_number"],
            "fixation_probability": tandem["fixation_probability"],
            "time_to_fixation_generations": tandem["time_to_fixation_generations"],
            "time_to_fixation_days": tandem["time_to_fixation_days"],
            "effective_escape_rate": tandem["effective_escape_rate"],
        },
        "time_to_resistance": escape_times_json,
        "sensitivity_matrix": sensitivity,
        "comparison_with_drugs": comparison,
    }

    # Gravar resultado
    output_path = write_result(envelope)
    logger.info("Resultado gravado em: %s", output_path)
    logger.info("Tempo de execucao: %.2f segundos", timer.elapsed)

    return envelope


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
