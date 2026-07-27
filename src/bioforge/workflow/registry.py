"""Step registry — maps workflow step ids to Python callables."""
from __future__ import annotations

from typing import Any, Callable

from bioforge.core.logging import get_logger

logger = get_logger("workflow.registry")


class StepRegistry:
    """A small dict-of-callables so workflows can name steps declaratively."""

    _instance: "StepRegistry | None" = None

    def __init__(self) -> None:
        self._steps: dict[str, Callable[..., Any]] = {}
        self._schemas: dict[str, dict[str, Any]] = {}

    @classmethod
    def instance(cls) -> "StepRegistry":
        if cls._instance is None:
            cls._instance = cls()
        return cls._instance

    def register(
        self,
        name: str,
        fn: Callable[..., Any],
        *,
        inputs: list[str] | None = None,
        outputs: list[str] | None = None,
        description: str = "",
    ) -> None:
        if name in self._steps:
            logger.warning("overwriting step registration '%s'", name)
        self._steps[name] = fn
        self._schemas[name] = {
            "inputs": list(inputs or []),
            "outputs": list(outputs or []),
            "description": description,
        }
        logger.debug("registered step '%s'", name)

    def get(self, name: str) -> Callable[..., Any]:
        if name not in self._steps:
            raise KeyError(f"unknown step '{name}' (registered: {sorted(self._steps)})")
        return self._steps[name]

    def schema(self, name: str) -> dict[str, Any]:
        return self._schemas.get(name, {})

    def known(self) -> list[str]:
        return sorted(self._steps)


_DEFAULT = StepRegistry.instance()


def register(name: str, *, inputs: list[str] | None = None,
             outputs: list[str] | None = None, description: str = ""):
    """Decorator: register a callable as a workflow step.

    The default registry (:func:`StepRegistry.instance`) is used by the
    executor and the CLI; tests can construct a fresh registry and pass
    it explicitly via :prop:`WorkflowExecutor.registry`.
    """

    def deco(fn: Callable[..., Any]) -> Callable[..., Any]:
        _DEFAULT.register(name, fn, inputs=inputs, outputs=outputs, description=description)
        return fn

    return deco
