"""Armazenamento e busca de embeddings vetoriais.

Fornece a classe EmbeddingStore para salvar, carregar e comparar
embeddings de proteinas, RNA e documentos. Usado pelos modulos:
    - 01_rag: embeddings de documentos para busca semantica
    - 03_leish_esm: embeddings proteicos do ESM-2
    - 04_rna_fm: embeddings de RNA
    - 07_contrastive: embeddings contrastivos epitopo-MHC

Formato de armazenamento: numpy .npy (eficiente e portavel).
Metadados ficam em .json sidecar com mesmo nome.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


class EmbeddingStore:
    """Gerencia armazenamento e busca de embeddings vetoriais.

    Cada store corresponde a um diretorio contendo pares de arquivos:
        - {name}.npy: matriz de embeddings (n_items x dim)
        - {name}.meta.json: metadados (ids, labels, parametros)

    Suporta busca por similaridade cosseno — suficiente para prototipo.
    Para producao, considerar FAISS ou Annoy.
    """

    def __init__(self, store_dir: Path) -> None:
        """Inicializa o store.

        Args:
            store_dir: Diretorio para armazenar embeddings.
        """
        self.store_dir = Path(store_dir)
        self.store_dir.mkdir(parents=True, exist_ok=True)

    def save(
        self,
        name: str,
        embeddings: np.ndarray,
        metadata: dict[str, Any] | None = None,
    ) -> Path:
        """Salva embeddings e metadados.

        Args:
            name: Nome do conjunto de embeddings (ex: "esm2_leish_proteome").
            embeddings: Matriz numpy (n_items x embedding_dim).
            metadata: Dict com informacoes adicionais (ids, labels, etc.).

        Returns:
            Caminho do arquivo .npy salvo.
        """
        npy_path = self.store_dir / f"{name}.npy"
        meta_path = self.store_dir / f"{name}.meta.json"

        np.save(npy_path, embeddings)

        meta = metadata or {}
        meta["shape"] = list(embeddings.shape)
        meta["dtype"] = str(embeddings.dtype)
        with open(meta_path, "w", encoding="utf-8") as fh:
            json.dump(meta, fh, indent=2, ensure_ascii=False)

        return npy_path

    def load(self, name: str) -> tuple[np.ndarray, dict[str, Any]]:
        """Carrega embeddings e metadados.

        Args:
            name: Nome do conjunto de embeddings.

        Returns:
            Tupla (embeddings, metadata).

        Raises:
            FileNotFoundError: Se o arquivo nao existir.
        """
        npy_path = self.store_dir / f"{name}.npy"
        meta_path = self.store_dir / f"{name}.meta.json"

        if not npy_path.exists():
            raise FileNotFoundError(
                f"Embeddings '{name}' nao encontrados em {npy_path}"
            )

        embeddings = np.load(npy_path)
        metadata: dict[str, Any] = {}
        if meta_path.exists():
            with open(meta_path, encoding="utf-8") as fh:
                metadata = json.load(fh)

        return embeddings, metadata

    def exists(self, name: str) -> bool:
        """Verifica se um conjunto de embeddings existe."""
        return (self.store_dir / f"{name}.npy").exists()

    def list_stores(self) -> list[str]:
        """Lista nomes de embeddings disponíveis neste diretorio."""
        return sorted(
            p.stem for p in self.store_dir.glob("*.npy")
        )

    @staticmethod
    def cosine_similarity(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Calcula similaridade cosseno entre vetores.

        Funciona para:
            - Vetor vs vetor: retorna escalar
            - Vetor vs matriz: retorna vetor de similaridades
            - Matriz vs matriz: retorna matriz de similaridades

        A similaridade cosseno mede o angulo entre vetores no espaco
        de embeddings. Valores proximos de 1.0 indicam alta similaridade.

        Args:
            a: Vetor ou matriz numpy.
            b: Vetor ou matriz numpy.

        Returns:
            Similaridade(s) cosseno no intervalo [-1, 1].
        """
        # Normaliza para vetores unitarios
        a_norm = a / (np.linalg.norm(a, axis=-1, keepdims=True) + 1e-10)
        b_norm = b / (np.linalg.norm(b, axis=-1, keepdims=True) + 1e-10)

        if a_norm.ndim == 1 and b_norm.ndim == 1:
            return np.dot(a_norm, b_norm)
        return a_norm @ b_norm.T
