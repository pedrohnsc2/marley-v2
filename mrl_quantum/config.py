"""Configuracao centralizada do pacote MRL-QUANTUM.

Importa parametros biologicos e termodinamicos do SSOT (aso_math.config)
e define constantes especificas dos modulos quanticos QAOA e VQE.

Hierarquia de dependencia:
    aso_math.config (SSOT) -> mrl_quantum.config -> mrl_quantum.qaoa / vqe
"""

from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

# ---------------------------------------------------------------------------
# Importacoes do SSOT — nunca duplicar estas constantes
# ---------------------------------------------------------------------------

from aso_math.config import ASO_SEQUENCE, SL_SEQUENCE, NN_PARAMS  # noqa: F401

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logger = logging.getLogger("mrl_quantum")

# ---------------------------------------------------------------------------
# Caminhos do projeto
# ---------------------------------------------------------------------------

PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parent.parent
MRL_QUANTUM_ROOT: Final[Path] = Path(__file__).resolve().parent
RESULTS_DIR: Final[Path] = MRL_QUANTUM_ROOT / "results"
QAOA_RESULTS_DIR: Final[Path] = RESULTS_DIR / "qaoa"
VQE_RESULTS_DIR: Final[Path] = RESULTS_DIR / "vqe"

# Caminhos dos certificados validados
MATH_CERT_PATH: Final[Path] = (
    PROJECT_ROOT / "aso_math" / "reports" / "results" / "math_certificate_v2.json"
)
DELIVERY_CERT_PATH: Final[Path] = (
    PROJECT_ROOT / "aso_delivery" / "reports" / "results" / "delivery_certificate.json"
)

# ---------------------------------------------------------------------------
# dG_BASE — carregado do certificado matematico (dimensao 1: wt_dg_kcal)
# Fallback: -27.97 kcal/mol (valor validado pelo pipeline)
# ---------------------------------------------------------------------------

_DG_FALLBACK: Final[float] = -27.97
_MATH_SCORE_FALLBACK: Final[str] = "52/60"
_DELIVERY_SCORE_FALLBACK: Final[str] = "60/60"


def _load_dg_base() -> tuple[float, str, str]:
    """Carrega dG_BASE do certificado matematico e scores dos certificados.

    Returns:
        Tupla (dg_base, math_score, delivery_score).
    """
    dg_base = _DG_FALLBACK
    math_score = _MATH_SCORE_FALLBACK
    delivery_score = _DELIVERY_SCORE_FALLBACK

    # Carregar dG do certificado matematico
    try:
        with open(MATH_CERT_PATH, encoding="utf-8") as fh:
            cert: dict[str, Any] = json.load(fh)
        # Dimensao 1 = THERMODYNAMIC OPTIMALITY
        for dim in cert.get("dimension_assessments", []):
            details = dim.get("details", {})
            if "wt_dg_kcal" in details:
                dg_base = float(details["wt_dg_kcal"])
                break
        math_score = f"{cert['composite_score']}/{cert['max_score']}"
        logger.info(
            "dG_BASE carregado do certificado: %.2f kcal/mol", dg_base
        )
    except (FileNotFoundError, KeyError, IndexError, json.JSONDecodeError) as exc:
        warnings.warn(
            f"Certificado matematico nao encontrado ou invalido ({exc}). "
            f"Usando fallback dG_BASE = {_DG_FALLBACK} kcal/mol.",
            stacklevel=2,
        )

    # Carregar score do certificado de delivery
    try:
        with open(DELIVERY_CERT_PATH, encoding="utf-8") as fh:
            dcert = json.load(fh)
        delivery_score = f"{dcert['composite_score']}/{dcert['max_score']}"
    except (FileNotFoundError, KeyError, json.JSONDecodeError) as exc:
        warnings.warn(
            f"Certificado de delivery nao encontrado ou invalido ({exc}). "
            f"Usando fallback delivery_score = {_DELIVERY_SCORE_FALLBACK}.",
            stacklevel=2,
        )

    return dg_base, math_score, delivery_score


DG_BASE: float
MATH_SCORE: str
DELIVERY_SCORE: str
DG_BASE, MATH_SCORE, DELIVERY_SCORE = _load_dg_base()

# ---------------------------------------------------------------------------
# Semente global para reproducibilidade
# ---------------------------------------------------------------------------

SEED: Final[int] = 42

# ---------------------------------------------------------------------------
# Perfil de execucao — "full" (padrao) ou "light" (Mac-friendly)
# ---------------------------------------------------------------------------
# "light" reduz camadas, iteracoes e shots para rodar em minutos sem
# estourar CPU/RAM em hardware tipico (MacBook).
# Selecionar via CLI: python -m mrl_quantum.run --profile light
# Ou via env: MRL_QUANTUM_PROFILE=light

import os as _os

PROFILE: Final[str] = _os.environ.get("MRL_QUANTUM_PROFILE", "full")

_IS_LIGHT: bool = PROFILE == "light"

# ---------------------------------------------------------------------------
# Parametros QAOA de alto nivel
# ---------------------------------------------------------------------------

QAOA_P_LAYERS: list[int] = [1] if _IS_LIGHT else [1, 2, 3]
QAOA_OPTIMIZERS: list[str] = ["COBYLA"] if _IS_LIGHT else ["COBYLA", "SPSA"]
QAOA_MAX_ITER: int = 50 if _IS_LIGHT else 300
QAOA_SHOTS: int = 1024 if _IS_LIGHT else 4096

# ---------------------------------------------------------------------------
# Configuracao QAOA — otimizacao de estereoisomeros PS
# ---------------------------------------------------------------------------
# O backbone fosforotioato (PS) do MRL-ASO-001 cria um centro quiral em
# cada ligacao internucleotidica. Para 25 nt existem 24 centros quirais,
# mas a otimizacao usa 25 qubits (1 por nucleotideo) para mapear tambem
# a conformacao do acucar adjacente.
#
# Correcoes energeticas baseadas em:
#   - Sp vs Rp: Wan WB et al. (2014) Nucleic Acids Res 42(22):13456-68
#   - Acoplamento J: Stec WJ et al. (2010) J Am Chem Soc 132(40):14133
#   - Penalidade Rp consecutivo: Eckstein F (2014) Nucleic Acids Ther 24(6)
# ---------------------------------------------------------------------------

# Constantes de estereoquimica PS como modulo-level Finals (acesso rapido)
SP_CORRECTION: Final[float] = -0.2       # kcal/mol por posicao Sp
RP_CORRECTION: Final[float] = +0.1       # kcal/mol por posicao Rp
J_COUPLING: Final[float] = -0.05         # cooperatividade entre posicoes adjacentes
RP_CONSECUTIVE_PENALTY: Final[float] = +0.5  # penalidade por 3+ Rp consecutivos

# Numero de posicoes PS (1 por nucleotideo do ASO)
N_POSITIONS: Final[int] = len(ASO_SEQUENCE)


@dataclass(frozen=True)
class QAOAConfig:
    """Parametros para o modulo QAOA de otimizacao estereoquimica PS."""

    # Correcao energetica Sp (kcal/mol) — estabiliza duplex com RNase H
    sp_correction: float = SP_CORRECTION

    # Correcao energetica Rp (kcal/mol) — menor afinidade, melhor uptake
    rp_correction: float = RP_CORRECTION

    # Acoplamento J entre centros quirais vizinhos (kcal/mol)
    j_coupling: float = J_COUPLING

    # Penalidade por 3+ Rp consecutivos (desestabiliza backbone) (kcal/mol)
    rp_consecutive_penalty: float = RP_CONSECUTIVE_PENALTY

    # Numero de qubits (1 por nucleotideo do ASO)
    n_qubits: int = N_POSITIONS

    # Camadas QAOA (p) — default conservador
    default_layers: int = 2

    # Numero de shots para amostragem
    n_shots: int = QAOA_SHOTS

    # Otimizador classico
    optimizer: str = "COBYLA"

    # Iteracoes maximas do otimizador
    max_iter: int = QAOA_MAX_ITER


# ---------------------------------------------------------------------------
# Configuracao VQE — simulacao do sitio ativo da RNase H1
# ---------------------------------------------------------------------------
# O VQE simula o sitio ativo da RNase H1 humana (PDB: 2QKB) para
# calcular a energia de ligacao do duplex ASO:RNA no sitio catalico.
# Usa active space restrito para viabilidade em simulador classico.
#
# Parametros baseados em:
#   - RNase H1 humana: Nowotny M et al. (2005) Mol Cell 17(5):691-700
#   - Active space: Reiher M et al. (2017) PNAS 114(29):7555-7560
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VQEConfig:
    """Parametros para o modulo VQE de simulacao do sitio ativo."""

    # Conjunto de base para calculo eletronico
    basis_set: str = "sto-3g"

    # Eletrons ativos no active space (ions Mg2+ do sitio catalico)
    active_electrons: int = 6

    # Orbitais ativos
    active_orbitals: int = 6

    # PDB do sitio ativo da RNase H1 humana
    pdb_id: str = "2QKB"

    # Ansatz VQE
    ansatz: str = "UCCSD"

    # Otimizador classico
    optimizer: str = "L-BFGS-B"

    # Iteracoes maximas
    max_iter: int = 50 if _IS_LIGHT else 200

    # Numero de shots
    n_shots: int = 1024 if _IS_LIGHT else 8192

    # Tolerancia de convergencia (Hartree)
    convergence_tol: float = 1e-6


# ---------------------------------------------------------------------------
# Backends disponiveis
# ---------------------------------------------------------------------------

BACKEND_CONFIGS: Final[dict[str, dict[str, Any]]] = {
    "aer_statevector": {
        "method": "statevector",
        "noise": False,
        "description": "Simulador exato sem ruido",
    },
    "aer_noisy": {
        "method": "density_matrix",
        "noise": True,
        "depolarizing_rate": 0.01,
        "description": "Simulador com modelo de ruido despolarizante",
    },
    "ibm_cloud": {
        "method": "cloud",
        "noise": True,
        "description": "Hardware IBM Quantum via qiskit_ibm_runtime",
    },
}

BACKENDS: Final[dict[str, str]] = {
    "aer_statevector": "Simulador statevector local (exato, sem ruido)",
    "aer_noisy": "Simulador Aer com modelo de ruido despolarizante",
    "ibm_cloud": "Hardware IBM Quantum via Qiskit Runtime (requer token)",
}

DEFAULT_BACKEND: Final[str] = "aer_statevector"

# Ordem de preferencia para fallback automatico
BACKEND_FALLBACK_ORDER: Final[list[str]] = [
    "ibm_cloud",
    "aer_noisy",
    "aer_statevector",
]
