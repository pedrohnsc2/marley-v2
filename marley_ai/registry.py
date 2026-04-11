"""Registro de modulos AI com decorador para auto-registro.

Segue o padrao de registro do aso_math.organisms, adaptado para
modulos de ML. Cada modulo se registra ao ser importado, permitindo
descoberta dinamica e execucao pelo orquestrador.

Uso:
    from marley_ai.registry import register, get_module, list_modules

    @register("03_leish_esm")
    class LeishESM:
        def configure(self, config): ...
        def validate_inputs(self): ...
        def run(self): ...
        def get_dependencies(self): ...
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


# ---------------------------------------------------------------------------
# Protocolo que todo modulo AI deve seguir
# ---------------------------------------------------------------------------

@runtime_checkable
class AIModule(Protocol):
    """Contrato que todo modulo do track AI deve implementar.

    O protocolo garante interface uniforme para o orquestrador,
    independente da implementacao interna de cada modulo.
    """

    def configure(self, config: Any) -> None:
        """Recebe configuracao e prepara o modulo para execucao.

        Args:
            config: AIModuleConfig ou subclasse especifica do modulo.
        """
        ...

    def validate_inputs(self) -> dict[str, Any]:
        """Verifica se todas as dependencias e inputs estao disponiveis.

        Returns:
            Dict com 'valid' (bool) e 'missing' (lista de dependencias ausentes).
        """
        ...

    def run(self) -> dict[str, Any]:
        """Executa o modulo e retorna resultados.

        Returns:
            Dict com resultados no formato envelope do marley_ai.
        """
        ...

    def get_dependencies(self) -> list[str]:
        """Lista slugs dos modulos que devem rodar antes deste.

        Returns:
            Lista de module_slugs (ex: ["01_rag", "03_leish_esm"]).
            Lista vazia se nao houver dependencias.
        """
        ...


# ---------------------------------------------------------------------------
# Registro global de modulos
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, type] = {}


def register(slug: str):
    """Decorador para registrar uma classe como modulo AI.

    Verifica em tempo de registro se a classe implementa o protocolo AIModule,
    prevenindo erros silenciosos em tempo de execucao.

    Args:
        slug: Identificador unico do modulo (ex: "03_leish_esm").

    Returns:
        Decorador que registra a classe e a retorna sem modificacao.

    Raises:
        ValueError: Se o slug ja estiver registrado (previne duplicatas).
        TypeError: Se a classe nao implementar o protocolo AIModule.
    """
    def decorator(cls: type) -> type:
        if slug in _REGISTRY:
            raise ValueError(
                f"Modulo '{slug}' ja registrado como {_REGISTRY[slug].__name__}. "
                f"Tentativa de registrar {cls.__name__}."
            )
        # Verificacao estrutural do protocolo — garante que a classe
        # tem todos os metodos necessarios antes de registrar.
        required_methods = ["configure", "validate_inputs", "run", "get_dependencies"]
        missing = [m for m in required_methods if not callable(getattr(cls, m, None))]
        if missing:
            raise TypeError(
                f"Classe {cls.__name__} nao implementa metodos obrigatorios: {missing}. "
                f"Veja marley_ai.registry.AIModule para o contrato."
            )
        _REGISTRY[slug] = cls
        return cls
    return decorator


def get_module(slug: str) -> type:
    """Retorna a classe registrada para um slug.

    Args:
        slug: Identificador do modulo.

    Returns:
        Classe registrada.

    Raises:
        KeyError: Se o slug nao estiver registrado.
    """
    if slug not in _REGISTRY:
        available = ", ".join(sorted(_REGISTRY.keys())) or "(nenhum registrado)"
        raise KeyError(
            f"Modulo '{slug}' nao encontrado no registro. Disponiveis: {available}"
        )
    return _REGISTRY[slug]


def list_modules() -> list[dict[str, str]]:
    """Retorna lista de modulos registrados com metadados.

    Returns:
        Lista de dicts com 'slug', 'class_name', e 'module_path'.
    """
    result = []
    for slug, cls in sorted(_REGISTRY.items()):
        result.append({
            "slug": slug,
            "class_name": cls.__name__,
            "module_path": f"{cls.__module__}.{cls.__qualname__}",
        })
    return result
