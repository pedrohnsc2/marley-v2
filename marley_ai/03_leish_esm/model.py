"""Wrapper para o modelo ESM-2 de linguagem proteica.

Encapsula o carregamento do ESM-2 via fair-esm (preferencial) ou
HuggingFace (fallback), com suporte a MPS (Apple Silicon), CUDA e CPU.

O ESM-2 (Lin et al., 2023, Science) e um modelo de linguagem proteica
treinado em ~250M de sequencias do UniRef50. A representacao aprendida
captura propriedades estruturais (contatos 3D, estrutura secundaria)
e funcionais (dominios, motifs) sem necessidade de alinhamento.

Estrategias de pooling:
    - "mean": media sobre residuos (excluindo tokens especiais)
      — melhor para similaridade entre sequencias
    - "cls": token [CLS] (primeiro token)
      — captura representacao global da sequencia
    - "per_residue": embedding de cada residuo individual
      — util para analise posicao-especifica (epitopos, binding sites)

Referencia:
    Lin Z et al. (2023) Evolutionary-scale prediction of atomic-level
    protein structure with a language model. Science 379(6637):1123-1130
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import torch

from .config import ESMConfig

logger = logging.getLogger(__name__)


class ESMModel:
    """Wrapper para inferencia com ESM-2.

    Carrega o modelo via fair-esm (torch hub) como estrategia primaria,
    com fallback para HuggingFace transformers. Suporta extracao de
    embeddings per-residue e per-sequence com diferentes pooling.

    Atributos:
        config: configuracao do modulo
        model: modelo ESM-2 carregado
        alphabet: alfabeto/tokenizer do ESM-2
        batch_converter: funcao para converter sequencias em tensors
        device: dispositivo de computacao (mps/cuda/cpu)
        backend: "fair-esm" ou "huggingface"
    """

    def __init__(self, config: ESMConfig) -> None:
        self.config = config
        self.model: Any = None
        self.alphabet: Any = None
        self.batch_converter: Any = None
        self.tokenizer: Any = None  # usado apenas no backend HuggingFace
        self.device = torch.device(config.device)
        self.backend: str = ""
        self._embedding_dim: int = config.embedding_dim

    def load(self) -> None:
        """Carrega o modelo ESM-2. Tenta fair-esm primeiro, fallback HuggingFace.

        O modelo e salvo em cache em MODELS_DIR (~2.5GB para o t33).
        Na primeira execucao, o download pode levar alguns minutos.

        Raises:
            RuntimeError: se nenhum backend conseguir carregar o modelo.
        """
        # Configura diretorio de cache para torch hub
        hub_dir = self.config.models_dir / "torch_hub"
        hub_dir.mkdir(parents=True, exist_ok=True)
        torch.hub.set_dir(str(hub_dir))

        try:
            self._load_fair_esm()
        except Exception as e:
            logger.warning(
                "fair-esm nao disponivel (%s). Tentando HuggingFace...", e
            )
            try:
                self._load_huggingface()
            except Exception as e2:
                raise RuntimeError(
                    f"Nenhum backend conseguiu carregar o ESM-2. "
                    f"fair-esm: {e} | HuggingFace: {e2}"
                ) from e2

        logger.info(
            "ESM-2 carregado via %s no device %s (dim=%d)",
            self.backend, self.device, self._embedding_dim,
        )

    def _load_fair_esm(self) -> None:
        """Carrega via biblioteca fair-esm (Facebook AI Research).

        Usa esm.pretrained para download e cache automatico.
        """
        import esm

        model_loader = getattr(esm.pretrained, self.config.esm_model)
        model, alphabet = model_loader()

        model = model.eval()
        # MPS pode ter problemas com alguns dtypes — float32 e mais seguro
        model = model.to(self.device)

        self.model = model
        self.alphabet = alphabet
        self.batch_converter = alphabet.get_batch_converter()
        self.backend = "fair-esm"

    def _load_huggingface(self) -> None:
        """Carrega via HuggingFace transformers (fallback).

        Usa AutoModel e AutoTokenizer com cache em models_dir.
        """
        from transformers import AutoModel, AutoTokenizer

        cache = str(self.config.models_dir / "huggingface")

        self.tokenizer = AutoTokenizer.from_pretrained(
            self.config.model_name,
            cache_dir=cache,
        )
        model = AutoModel.from_pretrained(
            self.config.model_name,
            cache_dir=cache,
        )
        model = model.eval().to(self.device)

        self.model = model
        self.backend = "huggingface"

    @torch.no_grad()
    def embed_sequences(
        self,
        sequences: list[tuple[str, str]],
        pool_strategy: str | None = None,
    ) -> dict[str, np.ndarray]:
        """Gera embeddings para uma lista de sequencias proteicas.

        Processa em batches para evitar OOM. Para cada sequencia, extrai
        o embedding da ultima camada do transformer e aplica pooling.

        Args:
            sequences: lista de (nome, sequencia_aminoacidica)
            pool_strategy: estrategia de pooling (override do config)

        Returns:
            Dict mapeando nome -> numpy array com embedding.
            Shape: (embedding_dim,) para mean/cls, (seq_len, embedding_dim)
            para per_residue.
        """
        if self.model is None:
            raise RuntimeError("Modelo nao carregado. Chame load() primeiro.")

        strategy = pool_strategy or self.config.pool_strategy
        results: dict[str, np.ndarray] = {}
        batch_size = self.config.batch_size

        # Truncar sequencias longas
        truncated = []
        for name, seq in sequences:
            if len(seq) > self.config.max_seq_length:
                logger.warning(
                    "Sequencia '%s' truncada de %d para %d residuos",
                    name, len(seq), self.config.max_seq_length,
                )
                seq = seq[:self.config.max_seq_length]
            truncated.append((name, seq))

        # Processar em batches
        for i in range(0, len(truncated), batch_size):
            batch = truncated[i:i + batch_size]

            if self.backend == "fair-esm":
                batch_embeddings = self._embed_fair_esm(batch, strategy)
            else:
                batch_embeddings = self._embed_huggingface(batch, strategy)

            results.update(batch_embeddings)

        return results

    def _embed_fair_esm(
        self,
        batch: list[tuple[str, str]],
        strategy: str,
    ) -> dict[str, np.ndarray]:
        """Extrai embeddings usando fair-esm.

        O batch_converter adiciona tokens especiais (<cls>, <eos>).
        Precisamos remove-los antes do pooling.
        """
        # batch_converter espera lista de (nome, sequencia)
        labels, strs, tokens = self.batch_converter(batch)
        tokens = tokens.to(self.device)

        # Forward pass — extrair representacoes de todas as camadas
        out = self.model(tokens, repr_layers=[self.config.layer])
        representations = out["representations"][self.config.layer]
        # Shape: (batch, seq_len_com_tokens_especiais, embedding_dim)

        results: dict[str, np.ndarray] = {}

        for idx, (name, seq) in enumerate(batch):
            # Remover tokens especiais: [CLS] no inicio e [EOS] no final
            # tokens especiais: posicao 0 = <cls>, posicao len+1 = <eos>
            seq_len = len(seq)
            residue_emb = representations[idx, 1:seq_len + 1, :]
            # Shape: (seq_len, embedding_dim)

            if strategy == "mean":
                emb = residue_emb.mean(dim=0).cpu().numpy()
            elif strategy == "cls":
                emb = representations[idx, 0, :].cpu().numpy()
            elif strategy == "per_residue":
                emb = residue_emb.cpu().numpy()
            else:
                raise ValueError(f"Estrategia de pooling invalida: '{strategy}'")

            results[name] = emb

        return results

    def _embed_huggingface(
        self,
        batch: list[tuple[str, str]],
        strategy: str,
    ) -> dict[str, np.ndarray]:
        """Extrai embeddings usando HuggingFace transformers.

        O tokenizer do HuggingFace trata cada aminoacido como um token
        separado (character-level tokenization para proteinas).
        """
        names = [name for name, _ in batch]
        seqs = [seq for _, seq in batch]

        # Tokenizar — cada aminoacido vira um token + tokens especiais
        encoded = self.tokenizer(
            seqs,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=self.config.max_seq_length + 2,  # +2 para tokens especiais
        )
        encoded = {k: v.to(self.device) for k, v in encoded.items()}

        # Forward pass
        out = self.model(**encoded)
        hidden = out.last_hidden_state
        # Shape: (batch, seq_len_com_padding_e_tokens_especiais, embedding_dim)

        attention_mask = encoded["attention_mask"]
        results: dict[str, np.ndarray] = {}

        for idx, name in enumerate(names):
            seq_len = len(seqs[idx])
            # Remover token [CLS] (pos 0) e [EOS] (ultima posicao real)
            residue_emb = hidden[idx, 1:seq_len + 1, :]

            if strategy == "mean":
                # Mean pooling considerando apenas tokens reais (sem padding)
                mask = attention_mask[idx, 1:seq_len + 1].unsqueeze(-1).float()
                emb = (residue_emb * mask).sum(dim=0) / mask.sum(dim=0).clamp(min=1)
                emb = emb.cpu().numpy()
            elif strategy == "cls":
                emb = hidden[idx, 0, :].cpu().numpy()
            elif strategy == "per_residue":
                emb = residue_emb.cpu().numpy()
            else:
                raise ValueError(f"Estrategia de pooling invalida: '{strategy}'")

            results[name] = emb

        return results

    def get_embedding_dim(self) -> int:
        """Retorna dimensao do embedding (1280 para t33, 480 para t12)."""
        return self._embedding_dim
