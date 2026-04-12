"""CLI entry point para o pacote MRL-QUANTUM.

Uso:
    python -m mrl_quantum.run                          # ambos modulos (QAOA + VQE)
    python -m mrl_quantum.run --module qaoa            # apenas QAOA
    python -m mrl_quantum.run --module vqe             # apenas VQE
    python -m mrl_quantum.run --module benchmark       # benchmark classico vs quantico
    python -m mrl_quantum.run --backend aer_noisy      # com modelo de ruido
    python -m mrl_quantum.run --layers 4               # QAOA com 4 camadas
    python -m mrl_quantum.run --basis 631g             # VQE com base 6-31G
    python -m mrl_quantum.run --backend ibm_cloud --token <TOKEN>

Cada modulo e independente; este script apenas orquestra a execucao
e delega a geracao de relatorios ao subpacote reporting.
"""

from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from typing import Any

from core.logger import get_logger

from mrl_quantum.config import (
    BACKENDS,
    DEFAULT_BACKEND,
    DG_BASE,
    DELIVERY_SCORE,
    MATH_SCORE,
    PROFILE,
    QAOAConfig,
    QAOA_MAX_ITER,
    QAOA_SHOTS,
    SEED,
    VQEConfig,
)

logger = get_logger("mrl_quantum")

# ---------------------------------------------------------------------------
# Mapeamento de nomes CLI -> basis_set interno
# ---------------------------------------------------------------------------

_BASIS_MAP: dict[str, str] = {
    "sto3g": "sto-3g",
    "sto-3g": "sto-3g",
    "631g": "6-31g",
    "6-31g": "6-31g",
}


def _build_parser() -> argparse.ArgumentParser:
    """Constroi o parser de argumentos CLI."""
    parser = argparse.ArgumentParser(
        prog="mrl_quantum",
        description=(
            "MRL-QUANTUM — Extensao quantica da pipeline Marley. "
            "Otimizacao estereoquimica (QAOA) e simulacao do sitio ativo (VQE)."
        ),
    )
    parser.add_argument(
        "--module",
        type=str,
        choices=["qaoa", "vqe", "both", "benchmark"],
        default="both",
        help="Modulo a executar (default: both).",
    )
    parser.add_argument(
        "--backend",
        type=str,
        choices=list(BACKENDS.keys()),
        default=DEFAULT_BACKEND,
        help=f"Backend quantico (default: {DEFAULT_BACKEND}).",
    )
    parser.add_argument(
        "--layers",
        type=int,
        default=None,
        help="Numero de camadas QAOA (default: 2).",
    )
    parser.add_argument(
        "--basis",
        type=str,
        choices=list(_BASIS_MAP.keys()),
        default="sto3g",
        help="Conjunto de base para VQE (default: sto3g).",
    )
    parser.add_argument(
        "--token",
        type=str,
        default=None,
        help="Token IBM Quantum (obrigatorio para backend ibm_cloud).",
    )
    parser.add_argument(
        "--profile",
        type=str,
        choices=["full", "light"],
        default=None,
        help=(
            "Perfil de execucao: 'light' reduz iteracoes/shots para "
            "rodar sem travar o Mac. Default: usa MRL_QUANTUM_PROFILE env ou 'full'."
        ),
    )
    return parser


# ---------------------------------------------------------------------------
# Spinner visual para o terminal
# ---------------------------------------------------------------------------

_SPINNER_FRAMES = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]


class Spinner:
    """Spinner animado com timer decorrido para operacoes longas."""

    def __init__(self, label: str) -> None:
        self._label = label
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._start_time = 0.0

    def start(self) -> "Spinner":
        self._start_time = time.time()
        self._thread = threading.Thread(target=self._spin, daemon=True)
        self._thread.start()
        return self

    def _spin(self) -> None:
        idx = 0
        while not self._stop.is_set():
            elapsed = time.time() - self._start_time
            mins, secs = divmod(int(elapsed), 60)
            frame = _SPINNER_FRAMES[idx % len(_SPINNER_FRAMES)]
            sys.stderr.write(
                f"\r  {frame} {self._label} [{mins:02d}:{secs:02d}]  "
            )
            sys.stderr.flush()
            idx += 1
            self._stop.wait(0.12)

    def stop(self, status: str = "ok") -> None:
        self._stop.set()
        if self._thread:
            self._thread.join()
        elapsed = time.time() - self._start_time
        mins, secs = divmod(int(elapsed), 60)
        icon = "✓" if status == "ok" else "✗"
        sys.stderr.write(
            f"\r  {icon} {self._label} [{mins:02d}:{secs:02d}]\n"
        )
        sys.stderr.flush()


def _validate_args(args: argparse.Namespace) -> None:
    """Valida combinacoes de argumentos."""
    if args.backend == "ibm_cloud" and not args.token:
        logger.error(
            "Backend ibm_cloud requer --token. "
            "Obtenha em https://quantum.ibm.com/account"
        )
        sys.exit(1)


def _run_qaoa(
    qaoa_config: QAOAConfig,
    backend: str,
    seed: int,
    is_light: bool = False,
) -> dict[str, Any] | None:
    """Executa o modulo QAOA de otimizacao estereoquimica."""
    logger.info("=" * 60)
    logger.info("Executando QAOA — Otimizacao de estereoisomeros PS")
    logger.info(
        "  qubits=%d, camadas=%d, max_iter=%d, shots=%d, backend=%s",
        qaoa_config.n_qubits,
        qaoa_config.default_layers,
        qaoa_config.max_iter,
        qaoa_config.n_shots,
        backend,
    )
    logger.info("=" * 60)

    spinner = Spinner("QAOA otimizando estereoisomeros PS...").start()
    start = time.time()
    try:
        from mrl_quantum.qaoa.qaoa_optimizer import run_qaoa  # type: ignore[import-not-found]

        result = run_qaoa(
            n_positions=qaoa_config.n_qubits,
            n_layers=qaoa_config.default_layers,
            optimizer_name=qaoa_config.optimizer if is_light else None,
            backend_name=backend,
            max_iter=qaoa_config.max_iter,
            seed=seed,
            max_qubits=10 if is_light else None,
        )
        elapsed = round(time.time() - start, 2)
        spinner.stop("ok")
        logger.info("QAOA concluido em %.2f s", elapsed)
        return result
    except ImportError:
        elapsed = round(time.time() - start, 2)
        spinner.stop("fail")
        logger.warning(
            "Modulo mrl_quantum.qaoa ainda nao implementado. "
            "Retornando None (%.2f s).",
            elapsed,
        )
        return None
    except Exception as exc:
        elapsed = round(time.time() - start, 2)
        spinner.stop("fail")
        logger.error("QAOA falhou (%.2f s): %s", elapsed, exc)
        return None


def _run_vqe(
    vqe_config: VQEConfig,
    backend: str,
    seed: int,
) -> dict[str, Any] | None:
    """Executa o modulo VQE de simulacao do sitio ativo."""
    logger.info("=" * 60)
    logger.info("Executando VQE — Simulacao do sitio ativo RNase H1")
    logger.info(
        "  base=%s, eletrons=%d, orbitais=%d, backend=%s",
        vqe_config.basis_set,
        vqe_config.active_electrons,
        vqe_config.active_orbitals,
        backend,
    )
    logger.info("=" * 60)

    spinner = Spinner("VQE simulando sitio ativo RNase H1...").start()
    start = time.time()
    try:
        from mrl_quantum.vqe.electronic_structure import run_vqe  # type: ignore[import-not-found]

        result = run_vqe(basis=vqe_config.basis_set, max_vqe_iter=vqe_config.max_iter)
        elapsed = round(time.time() - start, 2)
        spinner.stop("ok")
        logger.info("VQE concluido em %.2f s", elapsed)
        return result
    except ImportError:
        elapsed = round(time.time() - start, 2)
        spinner.stop("fail")
        logger.warning(
            "Modulo mrl_quantum.vqe ainda nao implementado. "
            "Retornando None (%.2f s).",
            elapsed,
        )
        return None
    except Exception as exc:
        elapsed = round(time.time() - start, 2)
        spinner.stop("fail")
        logger.error("VQE falhou (%.2f s): %s", elapsed, exc)
        return None


def _run_benchmark() -> None:
    """Executa benchmark comparativo classico vs quantico."""
    logger.info("=" * 60)
    logger.info("Benchmark classico vs quantico — ainda nao implementado")
    logger.info("=" * 60)
    logger.warning(
        "O modulo de benchmark sera implementado apos QAOA e VQE estarem prontos."
    )


def main(argv: list[str] | None = None) -> None:
    """Entry point do orquestrador MRL-QUANTUM."""
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)

    # Determinar profile efetivo: CLI > env > default ("full")
    profile = args.profile or PROFILE

    # --- Log de inicializacao ---
    logger.info("=" * 60)
    logger.info("MRL-QUANTUM v0.1.0 — Pipeline quantica do Marley")
    logger.info("=" * 60)
    logger.info("Profile: %s", profile)
    logger.info("Carregando parametros MRL-ASO-001 de aso_math/config.py")
    logger.info(
        "dG_BASE = %.2f kcal/mol (certificado math %s, delivery %s)",
        DG_BASE,
        MATH_SCORE,
        DELIVERY_SCORE,
    )
    logger.info("seed = %d", SEED)
    logger.info("backend = %s", args.backend)
    logger.info("")

    # --- Construir configs com overrides do CLI e profile ---
    is_light = profile == "light"

    qaoa_layers = args.layers if args.layers is not None else (1 if is_light else QAOAConfig().default_layers)
    qaoa_max_iter = 50 if is_light else QAOA_MAX_ITER
    qaoa_shots = 1024 if is_light else QAOA_SHOTS
    qaoa_config = QAOAConfig(
        default_layers=qaoa_layers,
        max_iter=qaoa_max_iter,
        n_shots=qaoa_shots,
    )

    basis = _BASIS_MAP.get(args.basis, "sto-3g")
    vqe_max_iter = 50 if is_light else VQEConfig().max_iter
    vqe_config = VQEConfig(basis_set=basis, max_iter=vqe_max_iter)

    if is_light:
        logger.info(
            "Modo LIGHT ativo: QAOA(layers=%d, iter=%d, shots=%d), "
            "VQE(iter=%d)",
            qaoa_layers, qaoa_max_iter, qaoa_shots, vqe_max_iter,
        )

    # --- Dispatch ---
    qaoa_result: dict[str, Any] | None = None
    vqe_result: dict[str, Any] | None = None

    if args.module == "benchmark":
        _run_benchmark()
        logger.info("Pipeline quantica concluida.")
        return

    if args.module in ("qaoa", "both"):
        qaoa_result = _run_qaoa(qaoa_config, args.backend, SEED, is_light=is_light)

    if args.module in ("vqe", "both"):
        vqe_result = _run_vqe(vqe_config, args.backend, SEED)

    # --- Geracao de relatorios ---
    if qaoa_result is not None or vqe_result is not None:
        logger.info("")
        spinner = Spinner("Gerando relatorios...").start()
        try:
            from mrl_quantum.reporting.generate_reports import generate_all_reports

            paths = generate_all_reports(
                qaoa_result=qaoa_result,
                vqe_result=vqe_result,
                backend=args.backend,
                seed=SEED,
            )
            spinner.stop("ok")
            for p in paths:
                logger.info("  -> %s", p)
        except Exception as exc:
            spinner.stop("fail")
            logger.error("Falha na geracao de relatorios: %s", exc)
    else:
        logger.info("")
        logger.info(
            "Nenhum resultado disponivel — relatorios nao gerados. "
            "Verifique se os modulos QAOA/VQE estao implementados."
        )

    logger.info("")
    logger.info("Pipeline quantica concluida.")


if __name__ == "__main__":
    main()
