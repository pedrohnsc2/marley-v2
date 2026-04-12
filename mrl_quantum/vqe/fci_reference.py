"""Referencia classica FCI para validacao do VQE.

Full Configuration Interaction (FCI) fornece a solucao exata para
o Hamiltoniano eletronico dentro de uma base finita. Serve como
benchmark para avaliar a acuracia do VQE.

Para sistemas pequenos (< 20 orbitais), FCI e computacionalmente
viavel e fornece a resposta "correta" contra a qual medimos o
erro do VQE.
"""

from __future__ import annotations

import logging
import time
from typing import Any, Final

import numpy as np

logger = logging.getLogger(__name__)

HARTREE_TO_KCAL: Final[float] = 627.509474
SEED: Final[int] = 42


def run_fci_reference(
    atom_string: str | None = None,
    basis: str = "sto3g",
    charge: int = 0,
    spin: int = 0,
    num_active_electrons: int = 6,
    num_active_orbitals: int = 6,
) -> dict[str, Any]:
    """Executa FCI classico via PySCF para o mesmo sistema do VQE.

    Pipeline:
    1. RHF para orbitais iniciais
    2. CASCI (FCI dentro do espaco ativo) para energia exata
    3. FCI completo (se viavel) para comparacao

    Args:
        atom_string: Geometria PySCF. Se None, usa modelo reduzido.
        basis: Conjunto de base.
        charge: Carga total.
        spin: 2S.
        num_active_electrons: Eletrons no espaco ativo para CASCI.
        num_active_orbitals: Orbitais no espaco ativo para CASCI.

    Returns:
        Dicionario com energias FCI, CASCI e comparacao.
    """
    from pyscf import fci, gto, mcscf, scf

    np.random.seed(SEED)
    t0 = time.time()

    # Modelo reduzido padrao se nao fornecido
    if atom_string is None:
        from mrl_quantum.vqe.electronic_structure import (
            _build_reduced_model_atom_string,
        )
        atom_string, charge, spin = _build_reduced_model_atom_string()
        logger.info("FCI usando modelo reduzido: %s", atom_string)

    # --- Passo 1: RHF ---
    mol = gto.M(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=0,
        symmetry=False,
    )
    mf = scf.RHF(mol)
    mf.max_cycle = 200
    mf.kernel()
    rhf_energy = mf.e_tot

    logger.info("FCI-ref RHF: %.8f Ha | convergiu: %s", rhf_energy, mf.converged)

    # --- Passo 2: CASCI (FCI no espaco ativo) ---
    total_e = mol.nelectron
    total_o = mf.mo_energy.shape[0]
    act_e = min(num_active_electrons, total_e)
    act_o = min(num_active_orbitals, total_o)

    # Ajustar para consistencia
    if act_e % 2 != 0:
        act_e -= 1
    if act_e > 2 * act_o:
        act_e = 2 * act_o

    cas = mcscf.CASCI(mf, act_o, act_e)
    e_casci = cas.kernel()[0]
    logger.info("CASCI(%de, %do): %.8f Ha", act_e, act_o, e_casci)

    # --- Passo 3: FCI completo (se sistema pequeno o suficiente) ---
    e_fci = None
    fci_feasible = total_o <= 16 and total_e <= 16

    if fci_feasible:
        try:
            cisolver = fci.FCI(mf)
            cisolver.max_cycle = 200
            cisolver.conv_tol = 1e-10
            e_fci_val, _ = cisolver.kernel()
            e_fci = float(e_fci_val)
            logger.info("FCI completo: %.8f Ha", e_fci)
        except Exception as exc:
            logger.warning("FCI completo falhou: %s", exc)
            e_fci = None
    else:
        logger.info(
            "FCI completo inviavel (%d orbitais, %d eletrons), usando CASCI",
            total_o,
            total_e,
        )

    elapsed = time.time() - t0

    # Melhor energia FCI disponivel
    best_fci = e_fci if e_fci is not None else e_casci

    result = {
        "rhf_energy_ha": round(float(rhf_energy), 8),
        "casci_energy_ha": round(float(e_casci), 8),
        "fci_energy_ha": round(float(e_fci), 8) if e_fci is not None else None,
        "best_fci_energy_ha": round(float(best_fci), 8),
        "best_fci_energy_kcal": round(float(best_fci * HARTREE_TO_KCAL), 4),
        "active_space": f"({act_e}e, {act_o}o)",
        "fci_feasible": fci_feasible,
        "correlation_energy_ha": round(float(best_fci - rhf_energy), 8),
        "correlation_energy_kcal": round(
            float((best_fci - rhf_energy) * HARTREE_TO_KCAL), 4
        ),
        "total_electrons": total_e,
        "total_orbitals": total_o,
        "converged": mf.converged,
        "runtime_seconds": round(elapsed, 2),
        "basis": basis,
        "atom_string": atom_string,
        "seed": SEED,
    }

    logger.info(
        "FCI referencia: melhor=%.8f Ha, correlacao=%.6f Ha (%.2f kcal/mol)",
        best_fci,
        best_fci - rhf_energy,
        (best_fci - rhf_energy) * HARTREE_TO_KCAL,
    )

    return result


def compare_vqe_fci(
    vqe_results: dict[str, Any],
    fci_results: dict[str, Any],
) -> dict[str, Any]:
    """Compara resultados VQE com referencia FCI.

    Args:
        vqe_results: Dicionario retornado por run_vqe().
        fci_results: Dicionario retornado por run_fci_reference().

    Returns:
        Dicionario com metricas de comparacao.
    """
    vqe_e = vqe_results["ground_state_energy_ha"]
    fci_e = fci_results["best_fci_energy_ha"]

    deviation_ha = vqe_e - fci_e
    deviation_mha = deviation_ha * 1000
    deviation_kcal = deviation_ha * HARTREE_TO_KCAL

    # Precisao quimica: < 1 kcal/mol (~1.6 mHa)
    chemical_accuracy = abs(deviation_kcal) < 1.0
    # Precisao computacional razoavel: < 10 mHa
    computational_accuracy = abs(deviation_mha) < 10.0

    comparison = {
        "vqe_energy_ha": vqe_e,
        "fci_energy_ha": fci_e,
        "deviation_ha": round(float(deviation_ha), 8),
        "deviation_mha": round(float(deviation_mha), 6),
        "deviation_kcal": round(float(deviation_kcal), 4),
        "chemical_accuracy": chemical_accuracy,
        "computational_accuracy": computational_accuracy,
        "vqe_method": vqe_results.get("method", "unknown"),
        "fci_active_space": fci_results["active_space"],
        "vqe_ansatz": vqe_results.get("ansatz", "unknown"),
        "vqe_num_qubits": vqe_results.get("num_qubits", 0),
    }

    logger.info(
        "VQE vs FCI: desvio=%.4f mHa (%.4f kcal/mol) | "
        "precisao quimica: %s | precisao computacional: %s",
        deviation_mha,
        deviation_kcal,
        chemical_accuracy,
        computational_accuracy,
    )

    return comparison
