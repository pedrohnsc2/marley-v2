"""VQE para estrutura eletronica do sitio ativo da RNase H1.

Pipeline completa:
1. PySCF: SCF (RHF) para orbitais iniciais
2. Espaco ativo: 6 eletrons em 6 orbitais (ou menor se necessario)
3. qiskit-nature: mapeamento para Hamiltoniano de qubits
4. VQE com ansatz UCCSD (ou EfficientSU2 como fallback)
5. Extracao de energia, gap HOMO-LUMO, cargas de Mulliken

Se a molecula do sitio ativo for grande demais para o mapeamento qubit
completo, usamos um modelo reduzido (fragmento Mg2+ + ligantes proximos)
que mantem a fisica essencial do mecanismo de dois-metais.

Referencia:
  Reiher M et al. (2017) PNAS 114(29):7555-7560 — analise de recursos
  quanticos para simulacao de nitrogenase (analogia com sitio bimetalico).
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any, Final

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------
HARTREE_TO_KCAL: Final[float] = 627.509474
SEED: Final[int] = 42
MAX_QUBITS: Final[int] = 20
DEFAULT_BASIS: Final[str] = "sto3g"

# Diretorio padrao para resultados
_DEFAULT_OUTPUT_DIR: Final[Path] = (
    Path(__file__).resolve().parent.parent.parent / "results" / "vqe"
)


def _atoms_to_pyscf_string(atoms: Any) -> str:
    """Converte ASE Atoms para string de geometria PySCF.

    Formato: 'Simbolo x y z; Simbolo x y z; ...'
    """
    parts = []
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        parts.append(f"{symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}")
    return "; ".join(parts)


def _compute_mulliken_charges(mf: Any) -> dict[str, list[float]]:
    """Calcula cargas de Mulliken a partir do objeto SCF convergido.

    Returns:
        Dicionario com simbolos atomicos e cargas correspondentes.
    """
    mol = mf.mol
    dm = mf.make_rdm1()
    s = mol.intor("int1e_ovlp")

    # Populacao de Mulliken: diag(D * S)
    pop = np.diag(dm @ s)

    # Cargas = Z_nucleo - populacao
    charges = []
    symbols = []
    for i in range(mol.natm):
        z = mol.atom_charge(i)
        sym = mol.atom_symbol(i)
        # Soma populacao dos orbitais deste atomo
        ao_slices = mol.aoslice_by_atom()
        start, end = ao_slices[i][2], ao_slices[i][3]
        atom_pop = np.sum(pop[start:end])
        charge = z - atom_pop
        charges.append(round(float(charge), 4))
        symbols.append(sym)

    return {"symbols": symbols, "charges": charges}


def _get_homo_lumo_gap(mf: Any) -> float:
    """Calcula gap HOMO-LUMO a partir dos autovalores orbitais (Ha).

    Para RHF, os orbitais estao em mf.mo_energy.
    """
    mo_occ = mf.mo_occ
    mo_energy = mf.mo_energy

    occupied = mo_energy[mo_occ > 0]
    virtual = mo_energy[mo_occ == 0]

    if len(occupied) == 0 or len(virtual) == 0:
        return 0.0

    homo = occupied[-1]
    lumo = virtual[0]
    return float(lumo - homo)


def _run_pyscf_scf(
    atom_string: str,
    basis: str = DEFAULT_BASIS,
    charge: int = 0,
    spin: int = 0,
) -> Any:
    """Executa RHF com PySCF.

    Args:
        atom_string: Geometria no formato PySCF.
        basis: Conjunto de base (padrao: STO-3G).
        charge: Carga total do sistema.
        spin: 2S (numero de eletrons desemparelhados).

    Returns:
        Objeto SCF convergido (mf).
    """
    from pyscf import gto, scf

    mol = gto.M(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=0,
        # Simetria desligada para sistemas grandes
        symmetry=False,
    )
    mf = scf.RHF(mol)
    mf.max_cycle = 200
    mf.conv_tol = 1e-8
    mf.kernel()

    if not mf.converged:
        logger.warning("RHF nao convergiu, tentando com DIIS mais agressivo")
        mf.diis_space = 12
        mf.max_cycle = 500
        mf.kernel()

    logger.info(
        "RHF energia: %.8f Ha | convergiu: %s | base: %s",
        mf.e_tot,
        mf.converged,
        basis,
    )
    return mf


def _build_reduced_model_atom_string() -> tuple[str, int, int]:
    """Constroi string atomica para modelo reduzido do sitio ativo.

    Modelo minimo que captura a quimica essencial:
    - 2 Mg2+ (ions catalíticos)
    - 2 H2O como ligantes (representam coordenacao)
    - 1 formato HCOO- (representa carboxilato D10 ponte)

    Carga total: +2 (2 Mg2+ = +4, formato = -1, sistema = +2 para ser viavel)
    Na verdade, usamos um modelo neutro e pequeno para viabilidade:

    Modelo: [Mg(H2O)]2+ simplificado como MgO com distancia de 2.1 A
    Para manter tratavel, usamos Mg2O2 como modelo bimetalico minimo.

    Returns:
        (atom_string, charge, spin)
    """
    # Modelo bimetalico minimo: Mg-O-Mg-O (anel de 4 membros)
    # Representa a ponte carboxilato entre os dois Mg2+
    # Geometria: anel planar com distancias experimentais
    mg_mg = 3.9
    mg_o = 2.1

    # Posicoes em geometria de losango
    atom_string = (
        f"Mg 0.0 0.0 0.0; "
        f"O {mg_o:.4f} 0.0 0.0; "
        f"Mg {mg_mg:.4f} 0.0 0.0; "
        f"O {mg_mg - mg_o:.4f} {mg_o * 0.8:.4f} 0.0"
    )
    # Carga 0, spin 0 — sistema neutro fechado
    return atom_string, 0, 0


def _run_vqe_qiskit_nature(
    atom_string: str,
    basis: str,
    charge: int,
    spin: int,
    num_active_electrons: int = 6,
    num_active_orbitals: int = 6,
    max_vqe_iter: int = 200,
) -> dict[str, Any]:
    """Executa VQE via pipeline qiskit-nature.

    Pipeline:
    1. PySCFDriver -> ElectronicStructureProblem
    2. ActiveSpaceTransformer (se necessario)
    3. ParityMapper com reducao de 2 qubits
    4. UCCSD ansatz + VQE

    Returns:
        Dicionario com resultados do VQE.
    """
    from qiskit.primitives import StatevectorEstimator
    from qiskit_algorithms import NumPyMinimumEigensolver, VQE
    from qiskit_algorithms.optimizers import COBYLA, L_BFGS_B
    from qiskit_nature.second_q.circuit.library import UCCSD, HartreeFock
    from qiskit_nature.second_q.drivers import MethodType, PySCFDriver
    from qiskit_nature.second_q.mappers import ParityMapper
    from qiskit_nature.second_q.transformers import ActiveSpaceTransformer

    t0 = time.time()
    convergence_history: list[float] = []

    # --- Passo 1: Driver PySCF ---
    driver = PySCFDriver(
        atom=atom_string,
        basis=basis,
        charge=charge,
        spin=spin,
        method=MethodType.RHF,
    )
    problem = driver.run()
    num_spatial = problem.num_spatial_orbitals
    num_particles = problem.num_particles
    nuclear_repulsion = problem.nuclear_repulsion_energy

    logger.info(
        "Problema eletronico: %d orbitais espaciais, particulas=%s, E_nuc=%.6f",
        num_spatial,
        num_particles,
        nuclear_repulsion,
    )

    # --- Passo 2: Espaco ativo (se necessario) ---
    total_electrons = sum(num_particles)
    need_active_space = (
        num_spatial > num_active_orbitals
        or total_electrons > num_active_electrons
    )

    if need_active_space:
        # Ajustar eletrons ativos para nao exceder os disponiveis
        actual_active_e = min(num_active_electrons, total_electrons)
        actual_active_o = min(num_active_orbitals, num_spatial)

        # Garantir que eletrons ativos sao par (RHF) e caibam nos orbitais
        if actual_active_e % 2 != 0:
            actual_active_e -= 1
        if actual_active_e > 2 * actual_active_o:
            actual_active_e = 2 * actual_active_o

        logger.info(
            "Aplicando ActiveSpaceTransformer: (%d e, %d o) de (%d e, %d o)",
            actual_active_e,
            actual_active_o,
            total_electrons,
            num_spatial,
        )

        ast = ActiveSpaceTransformer(
            num_electrons=actual_active_e,
            num_spatial_orbitals=actual_active_o,
        )
        problem = ast.transform(problem)
        num_spatial = problem.num_spatial_orbitals
        num_particles = problem.num_particles

    # --- Passo 3: Mapeamento para qubits ---
    # ParityMapper com reducao de 2 qubits (conserva numero de particulas)
    mapper = ParityMapper(num_particles=num_particles)

    fermionic_op = problem.hamiltonian.second_q_op()
    qubit_op = mapper.map(fermionic_op)
    num_qubits = qubit_op.num_qubits

    logger.info(
        "Hamiltoniano de qubits: %d qubits (limite: %d)",
        num_qubits,
        MAX_QUBITS,
    )

    if num_qubits > MAX_QUBITS:
        logger.warning(
            "Numero de qubits (%d) excede limite (%d). "
            "Reduzindo espaco ativo.",
            num_qubits,
            MAX_QUBITS,
        )
        # Reducao agressiva: 4 eletrons em 4 orbitais
        ast2 = ActiveSpaceTransformer(
            num_electrons=4, num_spatial_orbitals=4
        )
        # Re-rodar driver
        problem = driver.run()
        problem = ast2.transform(problem)
        num_spatial = problem.num_spatial_orbitals
        num_particles = problem.num_particles
        mapper = ParityMapper(num_particles=num_particles)
        fermionic_op = problem.hamiltonian.second_q_op()
        qubit_op = mapper.map(fermionic_op)
        num_qubits = qubit_op.num_qubits
        logger.info("Apos reducao: %d qubits", num_qubits)

    # --- Passo 4: Referencia classica (NumPy exato) ---
    numpy_solver = NumPyMinimumEigensolver()
    result_exact = numpy_solver.compute_minimum_eigenvalue(qubit_op)
    exact_energy = result_exact.eigenvalue.real
    logger.info("Energia exata (NumPy): %.8f Ha", exact_energy)

    # --- Passo 5: VQE ---
    # Ansatz: UCCSD (primeira opcao)
    ansatz_type = "UCCSD"
    try:
        hf_state = HartreeFock(num_spatial, num_particles, mapper)
        ansatz = UCCSD(
            num_spatial, num_particles, mapper, initial_state=hf_state
        )
        num_params = ansatz.num_parameters
        logger.info("Ansatz UCCSD: %d parametros", num_params)

        # Se muitos parametros, usar EfficientSU2
        if num_params > 500:
            raise ValueError(
                f"UCCSD tem {num_params} parametros, usando EfficientSU2"
            )

    except Exception as exc:
        logger.warning("UCCSD indisponivel (%s), usando EfficientSU2", exc)
        from qiskit.circuit.library import EfficientSU2

        ansatz = EfficientSU2(num_qubits, reps=2, entanglement="circular")
        ansatz_type = "EfficientSU2"
        num_params = ansatz.num_parameters

    # Callback para registrar convergencia
    def _callback(eval_count: int, parameters: np.ndarray, value: float,
                  metadata: dict | None) -> None:
        convergence_history.append(float(value))

    # Estimator e otimizador
    estimator = StatevectorEstimator(seed=SEED)
    initial_point = np.zeros(num_params)

    # COBYLA como primario (derivative-free, padrao em VQE)
    # L-BFGS-B como fallback (baseado em gradiente)
    optimizer_name = "COBYLA"
    try:
        optimizer = COBYLA(maxiter=max_vqe_iter)
        vqe = VQE(estimator, ansatz, optimizer, callback=_callback)
        vqe.initial_point = initial_point
        result_vqe = vqe.compute_minimum_eigenvalue(qubit_op)
    except Exception as exc:
        logger.warning("COBYLA falhou (%s), tentando L-BFGS-B", exc)
        optimizer_name = "L_BFGS_B"
        convergence_history.clear()
        optimizer = L_BFGS_B(maxiter=max_vqe_iter)
        vqe = VQE(estimator, ansatz, optimizer, callback=_callback)
        vqe.initial_point = initial_point
        result_vqe = vqe.compute_minimum_eigenvalue(qubit_op)

    vqe_eigenvalue = result_vqe.eigenvalue.real
    elapsed = time.time() - t0

    # --- Passo 6: Resultados ---
    # Usar problem.interpret() para obter energia total correta
    # (inclui energia inativa dos eletrons congelados pelo ActiveSpaceTransformer)
    vqe_interpreted = problem.interpret(result_vqe)
    total_vqe = float(vqe_interpreted.total_energies[0].real)

    exact_interpreted = problem.interpret(result_exact)
    total_exact = float(exact_interpreted.total_energies[0].real)

    # Desvio entre VQE e exato (no mesmo espaco ativo)
    deviation_ha = abs(total_vqe - total_exact)
    deviation_mha = deviation_ha * 1000

    logger.info(
        "VQE eigenvalue: %.8f Ha | total (com inativo): %.8f Ha | "
        "exato total: %.8f Ha | desvio: %.4f mHa | tempo: %.1f s",
        vqe_eigenvalue,
        total_vqe,
        total_exact,
        deviation_mha,
        elapsed,
    )

    return {
        "ground_state_energy_ha": round(float(total_vqe), 8),
        "ground_state_energy_kcal": round(float(total_vqe * HARTREE_TO_KCAL), 4),
        "electronic_energy_ha": round(float(vqe_eigenvalue), 8),
        "nuclear_repulsion_ha": round(float(nuclear_repulsion), 8),
        "exact_energy_ha": round(float(total_exact), 8),
        "deviation_mha": round(float(deviation_mha), 6),
        "deviation_percent": round(
            float(deviation_ha / abs(total_exact) * 100) if total_exact != 0 else 0.0,
            6,
        ),
        "num_qubits": int(num_qubits),
        "num_parameters": int(num_params),
        "ansatz": ansatz_type,
        "optimizer": optimizer_name,
        "num_particles": list(num_particles),
        "num_spatial_orbitals": int(num_spatial),
        "convergence_history": convergence_history,
        "converged": deviation_mha < 50.0,  # < 50 mHa = convergido (realista para UCCSD/COBYLA)
        "runtime_seconds": round(elapsed, 2),
        "seed": SEED,
    }


def _run_pyscf_only(
    atom_string: str,
    basis: str,
    charge: int,
    spin: int,
    num_active_electrons: int = 6,
    num_active_orbitals: int = 6,
) -> dict[str, Any]:
    """Fallback: executa SCF + CASCI com PySCF puro (sem qiskit).

    Usado quando o pipeline qiskit-nature falha.
    Ainda fornece resultados quimicos validos.
    """
    from pyscf import gto, mcscf, scf

    t0 = time.time()

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

    # Ajustar espaco ativo
    total_e = mol.nelectron
    total_o = mf.mo_energy.shape[0]
    act_e = min(num_active_electrons, total_e)
    act_o = min(num_active_orbitals, total_o)
    if act_e % 2 != 0:
        act_e -= 1
    if act_e > 2 * act_o:
        act_e = 2 * act_o

    # CASCI
    cas = mcscf.CASCI(mf, act_o, act_e)
    e_cas = cas.kernel()[0]

    # Gap e cargas
    gap = _get_homo_lumo_gap(mf)
    charges_data = _compute_mulliken_charges(mf)

    elapsed = time.time() - t0

    logger.info(
        "PySCF-only: RHF=%.8f Ha, CASCI=%.8f Ha, gap=%.4f Ha, tempo=%.1f s",
        mf.e_tot,
        e_cas,
        gap,
        elapsed,
    )

    return {
        "ground_state_energy_ha": round(float(e_cas), 8),
        "ground_state_energy_kcal": round(float(e_cas * HARTREE_TO_KCAL), 4),
        "electronic_energy_ha": round(float(e_cas - mol.energy_nuc()), 8),
        "nuclear_repulsion_ha": round(float(mol.energy_nuc()), 8),
        "rhf_energy_ha": round(float(mf.e_tot), 8),
        "exact_energy_ha": round(float(e_cas), 8),
        "deviation_mha": 0.0,  # CASCI e a referencia neste caso
        "deviation_percent": 0.0,
        "homo_lumo_gap_ha": round(float(gap), 6),
        "homo_lumo_gap_ev": round(float(gap * 27.2114), 4),
        "mulliken_charges": charges_data,
        "num_qubits": 0,
        "num_parameters": 0,
        "ansatz": "CASCI_pyscf",
        "optimizer": "direct_diag",
        "active_space": f"({act_e}e, {act_o}o)",
        "converged": mf.converged,
        "runtime_seconds": round(elapsed, 2),
        "method": "pyscf_casci",
        "seed": SEED,
    }


def run_vqe(
    atoms: Any | None = None,
    basis: str = DEFAULT_BASIS,
    charge: int | None = None,
    spin: int = 0,
    num_active_electrons: int = 4,
    num_active_orbitals: int = 4,
    max_vqe_iter: int = 500,
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Funcao principal: executa VQE para o sitio ativo da RNase H1.

    Pipeline com fallbacks:
    1. Tenta pipeline completa qiskit-nature (PySCFDriver + VQE)
    2. Se falhar, executa PySCF puro com CASCI

    Espaco ativo padrao: 4 eletrons em 4 orbitais (6 qubits apos ParityMapper).
    Para calculos mais acurados, usar num_active_electrons=6, num_active_orbitals=6
    (10 qubits, ~117 parametros UCCSD, runtime ~10min).

    Args:
        atoms: Objeto ASE Atoms. Se None, constroi modelo reduzido.
        basis: Conjunto de base (padrao: sto3g).
        charge: Carga total. Se None, detecta automaticamente.
        spin: 2S (padrao: 0 para sistema de camada fechada).
        num_active_electrons: Eletrons no espaco ativo (padrao: 4).
        num_active_orbitals: Orbitais no espaco ativo (padrao: 4).
        max_vqe_iter: Maximo de iteracoes VQE (padrao: 500 para COBYLA).
        output_dir: Diretorio para resultados.

    Returns:
        Dicionario com resultados completos do calculo.
    """
    np.random.seed(SEED)
    out = output_dir or _DEFAULT_OUTPUT_DIR
    out.mkdir(parents=True, exist_ok=True)

    # Se nao recebeu atoms, construir modelo reduzido
    if atoms is not None:
        atom_string = _atoms_to_pyscf_string(atoms)
        if charge is None:
            charge = 0
    else:
        atom_string, charge, spin = _build_reduced_model_atom_string()
        logger.info("Usando modelo reduzido: %s", atom_string)

    # --- Fase 1: SCF para diagnostico e cargas ---
    mf = _run_pyscf_scf(atom_string, basis, charge, spin)
    gap = _get_homo_lumo_gap(mf)
    charges_data = _compute_mulliken_charges(mf)

    logger.info("HOMO-LUMO gap: %.4f Ha (%.2f eV)", gap, gap * 27.2114)

    # --- Fase 2: Tentar pipeline VQE qiskit-nature ---
    results: dict[str, Any]
    try:
        results = _run_vqe_qiskit_nature(
            atom_string=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            num_active_electrons=num_active_electrons,
            num_active_orbitals=num_active_orbitals,
            max_vqe_iter=max_vqe_iter,
        )
        results["method"] = "vqe_qiskit_nature"
    except Exception as exc:
        logger.warning(
            "Pipeline qiskit-nature falhou (%s), usando PySCF puro", exc
        )
        results = _run_pyscf_only(
            atom_string=atom_string,
            basis=basis,
            charge=charge,
            spin=spin,
            num_active_electrons=num_active_electrons,
            num_active_orbitals=num_active_orbitals,
        )

    # Adicionar informacoes extras
    results["homo_lumo_gap_ha"] = round(float(gap), 6)
    results["homo_lumo_gap_ev"] = round(float(gap * 27.2114), 4)
    results["mulliken_charges"] = charges_data
    results["basis"] = basis
    results["atom_string"] = atom_string

    # Salvar convergencia como JSON
    conv_path = out / "vqe_convergence.json"
    with open(conv_path, "w") as fh:
        json.dump(
            {
                "convergence_history": results.get("convergence_history", []),
                "final_energy_ha": results["ground_state_energy_ha"],
                "method": results.get("method", "unknown"),
            },
            fh,
            indent=2,
        )

    logger.info(
        "VQE concluido: E=%.8f Ha (%.2f kcal/mol) | metodo: %s",
        results["ground_state_energy_ha"],
        results["ground_state_energy_kcal"],
        results.get("method", "unknown"),
    )

    return results
