"""Testes VQE usando molecula H2 (caso mais simples possivel).

H2 com distancia de ligacao 0.735 A na base STO-3G e o sistema
benchmark padrao para validacao de algoritmos VQE. A energia
FCI exata e -1.1373 Ha.

Todos os testes devem rodar em < 30s no CPU.
"""

from __future__ import annotations

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Constantes do benchmark H2
# ---------------------------------------------------------------------------
H2_BOND_LENGTH: float = 0.735  # Angstrom
H2_ATOM_STRING: str = f"H 0 0 0; H 0 0 {H2_BOND_LENGTH}"
H2_BASIS: str = "sto3g"
H2_FCI_ENERGY: float = -1.1373  # Ha (referencia literatura STO-3G)
H2_TOLERANCE: float = 0.05     # Ha — tolerancia para VQE


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def h2_fci_results() -> dict:
    """Executa FCI para H2 uma vez e compartilha entre testes."""
    from mrl_quantum.vqe.fci_reference import run_fci_reference

    return run_fci_reference(
        atom_string=H2_ATOM_STRING,
        basis=H2_BASIS,
        charge=0,
        spin=0,
        num_active_electrons=2,
        num_active_orbitals=2,
    )


@pytest.fixture(scope="module")
def h2_vqe_results() -> dict:
    """Executa VQE para H2 uma vez e compartilha entre testes."""
    from mrl_quantum.vqe.electronic_structure import run_vqe
    from ase import Atoms

    # Construir H2 com ASE
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, H2_BOND_LENGTH]])

    return run_vqe(
        atoms=h2,
        basis=H2_BASIS,
        charge=0,
        spin=0,
        num_active_electrons=2,
        num_active_orbitals=2,
        max_vqe_iter=100,
    )


# ---------------------------------------------------------------------------
# Testes
# ---------------------------------------------------------------------------
class TestH2FCI:
    """Testes para FCI classico com H2."""

    def test_fci_energy_matches_literature(
        self, h2_fci_results: dict
    ) -> None:
        """Energia FCI deve estar proxima do valor de literatura (-1.1373 Ha)."""
        best_fci = h2_fci_results["best_fci_energy_ha"]
        assert abs(best_fci - H2_FCI_ENERGY) < 0.001, (
            f"FCI energia {best_fci:.6f} Ha difere do esperado "
            f"{H2_FCI_ENERGY:.4f} Ha por mais de 1 mHa"
        )

    def test_fci_lower_than_rhf(self, h2_fci_results: dict) -> None:
        """Energia FCI deve ser menor que RHF (correlacao e negativa)."""
        assert h2_fci_results["best_fci_energy_ha"] < h2_fci_results["rhf_energy_ha"], (
            "Energia de correlacao deve ser negativa (FCI < RHF)"
        )

    def test_fci_correlation_energy(self, h2_fci_results: dict) -> None:
        """Energia de correlacao para H2/STO-3G deve ser pequena (~-20 mHa)."""
        corr = h2_fci_results["correlation_energy_ha"]
        # Correlacao negativa e pequena para H2
        assert corr < 0, "Correlacao deve ser negativa"
        assert abs(corr) < 0.1, (
            f"Correlacao {corr:.6f} Ha muito grande para H2/STO-3G"
        )


class TestH2VQE:
    """Testes para VQE com H2."""

    def test_vqe_energy_near_expected(self, h2_vqe_results: dict) -> None:
        """VQE deve encontrar energia proxima de -1.137 Ha (dentro de 0.05 Ha)."""
        vqe_e = h2_vqe_results["ground_state_energy_ha"]
        assert abs(vqe_e - H2_FCI_ENERGY) < H2_TOLERANCE, (
            f"VQE energia {vqe_e:.6f} Ha difere do esperado "
            f"{H2_FCI_ENERGY:.4f} Ha por mais de {H2_TOLERANCE} Ha"
        )

    def test_vqe_converged(self, h2_vqe_results: dict) -> None:
        """VQE deve convergir (desvio < 10 mHa do exato)."""
        assert h2_vqe_results["converged"], (
            f"VQE nao convergiu. Desvio: {h2_vqe_results['deviation_mha']:.4f} mHa"
        )

    def test_vqe_vs_fci_comparison(
        self, h2_vqe_results: dict, h2_fci_results: dict
    ) -> None:
        """VQE e FCI devem concordar dentro de 10 mHa para H2."""
        from mrl_quantum.vqe.fci_reference import compare_vqe_fci

        comparison = compare_vqe_fci(h2_vqe_results, h2_fci_results)
        assert comparison["computational_accuracy"], (
            f"VQE vs FCI desvio {comparison['deviation_mha']:.4f} mHa > 10 mHa"
        )

    def test_vqe_has_convergence_history(self, h2_vqe_results: dict) -> None:
        """VQE deve registrar historico de convergencia."""
        history = h2_vqe_results.get("convergence_history", [])
        assert len(history) > 0, "Historico de convergencia vazio"
        # Energia final deve ser menor ou igual a inicial
        # (L-BFGS-B pode nao ser monotono, mas deve melhorar no geral)
        assert history[-1] <= history[0] + 0.1, (
            "Energia final deve ser proxima ou menor que inicial"
        )

    def test_vqe_metadata_complete(self, h2_vqe_results: dict) -> None:
        """Resultados VQE devem conter todos os campos esperados."""
        required_keys = [
            "ground_state_energy_ha",
            "ground_state_energy_kcal",
            "num_qubits",
            "ansatz",
            "optimizer",
            "seed",
            "runtime_seconds",
            "homo_lumo_gap_ha",
            "mulliken_charges",
            "basis",
        ]
        for key in required_keys:
            assert key in h2_vqe_results, f"Campo '{key}' ausente nos resultados"

    def test_vqe_seed_determinism(self) -> None:
        """Duas execucoes com mesmo seed devem dar mesmo resultado."""
        from mrl_quantum.vqe.electronic_structure import run_vqe
        from ase import Atoms

        h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, H2_BOND_LENGTH]])

        r1 = run_vqe(
            atoms=h2, basis=H2_BASIS, charge=0, spin=0,
            num_active_electrons=2, num_active_orbitals=2, max_vqe_iter=50,
        )
        r2 = run_vqe(
            atoms=h2, basis=H2_BASIS, charge=0, spin=0,
            num_active_electrons=2, num_active_orbitals=2, max_vqe_iter=50,
        )

        assert abs(r1["ground_state_energy_ha"] - r2["ground_state_energy_ha"]) < 1e-6, (
            "Resultados com mesmo seed devem ser deterministicos"
        )
