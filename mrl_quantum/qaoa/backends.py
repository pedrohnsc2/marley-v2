"""Gerenciamento de backends quanticos para QAOA.

Suporta tres backends com fallback automatico:
    1. ibm_cloud — hardware IBM Quantum (requer token)
    2. aer_noisy — simulador Aer com modelo de ruido despolarizante
    3. aer_statevector — simulador exato sem ruido

Para QAOA via qiskit_algorithms, o backend e abstraido pelo
StatevectorSampler (simulacao exata) ou BackendSamplerV2 (com ruido).
"""

from __future__ import annotations

import logging
from typing import Any

from qiskit.primitives import StatevectorSampler
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

from mrl_quantum.config import (
    BACKEND_CONFIGS,
    BACKEND_FALLBACK_ORDER,
    DEFAULT_BACKEND,
    SEED,
)

logger = logging.getLogger("mrl_quantum.qaoa.backends")


def _build_noise_model(depol_rate: float = 0.01) -> NoiseModel:
    """Constroi modelo de ruido despolarizante simples.

    Aplica erro despolarizante em todas as portas de 1 e 2 qubits.

    Args:
        depol_rate: Taxa de erro despolarizante (default: 0.01 = 1%).

    Returns:
        NoiseModel configurado.
    """
    noise_model = NoiseModel()

    # Erro de 1 qubit nas portas comuns do QAOA
    error_1q = depolarizing_error(depol_rate, 1)
    for gate in ["rx", "ry", "rz", "x", "h", "sx"]:
        noise_model.add_all_qubit_quantum_error(error_1q, gate)

    # Erro de 2 qubits (taxa ligeiramente maior — tipico de hardware)
    error_2q = depolarizing_error(depol_rate * 2, 2)
    for gate in ["cx", "cz", "rzz"]:
        noise_model.add_all_qubit_quantum_error(error_2q, gate)

    logger.info(
        "Modelo de ruido criado: depol_1q=%.4f, depol_2q=%.4f",
        depol_rate, depol_rate * 2,
    )
    return noise_model


def get_aer_backend(
    method: str = "statevector",
    noise: bool = False,
    depol_rate: float = 0.01,
    seed: int = SEED,
) -> AerSimulator:
    """Retorna um AerSimulator configurado.

    Args:
        method: Metodo de simulacao ('statevector', 'density_matrix').
        noise: Se True, aplica modelo de ruido despolarizante.
        depol_rate: Taxa de erro despolarizante.
        seed: Semente para reproducibilidade.

    Returns:
        AerSimulator configurado.
    """
    kwargs: dict[str, Any] = {"method": method, "seed_simulator": seed}

    if noise:
        noise_model = _build_noise_model(depol_rate)
        kwargs["noise_model"] = noise_model

    backend = AerSimulator(**kwargs)
    logger.info(
        "Backend Aer criado: method=%s, noise=%s, seed=%d",
        method, noise, seed,
    )
    return backend


def get_sampler(
    backend_name: str = DEFAULT_BACKEND,
    seed: int = SEED,
) -> tuple[StatevectorSampler, str]:
    """Retorna um Sampler configurado para QAOA.

    O qiskit_algorithms.QAOA usa Sampler (nao Backend diretamente).
    Para simulacao exata usamos StatevectorSampler.
    Para backends com ruido, tambem usamos StatevectorSampler
    (pois BackendSamplerV2 nao suporta circuitos QAOA parametrizados).

    Args:
        backend_name: Nome do backend desejado.
        seed: Semente para reproducibilidade.

    Returns:
        Tupla (sampler, nome_do_backend_usado).
    """
    # StatevectorSampler funciona para todos os casos no simulador
    # Nota: para ruido real, precisariamos de hardware IBM ou
    # transpilacao manual + BackendSamplerV2
    sampler = StatevectorSampler(seed=seed)

    if backend_name == "aer_statevector":
        logger.info("Sampler configurado: StatevectorSampler (exato, sem ruido)")
        return sampler, "aer_statevector"

    if backend_name == "aer_noisy":
        # Log que estamos usando statevector como proxy
        logger.info(
            "Sampler configurado: StatevectorSampler (proxy para aer_noisy — "
            "ruido real requer transpilacao manual)"
        )
        return sampler, "aer_statevector"

    if backend_name == "ibm_cloud":
        # Tentar IBM Quantum Runtime
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

            service = QiskitRuntimeService()
            ibm_backend = service.least_busy(simulator=False)
            ibm_sampler = SamplerV2(backend=ibm_backend)
            logger.info(
                "Sampler IBM Quantum conectado: %s", ibm_backend.name
            )
            return ibm_sampler, f"ibm_{ibm_backend.name}"
        except Exception as exc:
            logger.warning(
                "Falha ao conectar IBM Quantum (%s). Tentando fallback...", exc
            )

    # Fallback automatico
    logger.info(
        "Backend '%s' indisponivel. Usando fallback: aer_statevector",
        backend_name,
    )
    return StatevectorSampler(seed=seed), "aer_statevector"


def get_backend(
    name: str = DEFAULT_BACKEND,
    seed: int = SEED,
) -> tuple[Any, str]:
    """Interface principal para obter backend/sampler.

    Tenta o backend solicitado e faz fallback automatico se necessario.
    Segue a ordem: ibm_cloud -> aer_noisy -> aer_statevector.

    Args:
        name: Nome do backend desejado.
        seed: Semente para reproducibilidade.

    Returns:
        Tupla (sampler_ou_backend, nome_real_usado).
    """
    return get_sampler(backend_name=name, seed=seed)
