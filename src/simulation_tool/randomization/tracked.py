import json
from functools import wraps
from pathlib import Path
from typing import Any, Callable, NewType, ParamSpec, TypedDict, TypeVar

P = ParamSpec("P")
R = TypeVar("R")

VariableName = NewType("VariableName", str)


class RandomizationEntry(TypedDict):
    method: str
    args: tuple
    kwargs: dict[str, Any]


RANDOMIZATION_LOG: dict[VariableName, RandomizationEntry] = {}


def record(
    variable_name: str,
    method: str,
    args: tuple,
    kwargs: dict[str, Any],
):
    if variable_name in RANDOMIZATION_LOG:
        return

    RANDOMIZATION_LOG[variable_name] = {
        "method": method,
        "args": args,
        "kwargs": kwargs,
    }


def reset():
    RANDOMIZATION_LOG.clear()


def export_json(path: Path):
    with open(path, "w") as f:
        json.dump(RANDOMIZATION_LOG, f, indent=2)


def tracked(*, is_method: bool = False):
    def decorator(func: Callable[P, R]):
        name = func.__qualname__

        if is_method:

            @wraps(func)
            def wrapper(
                self, variable: VariableName, *args: P.args, **kwargs: P.kwargs
            ) -> R:
                value = func(self, *args, **kwargs)
                record(variable, name, args, kwargs)
                return value

            return wrapper

        @wraps(func)
        def wrapper(variable: VariableName, *args: P.args, **kwargs: P.kwargs) -> R:
            value = func(*args, **kwargs)
            record(variable, name, args, kwargs)
            return value

        return wrapper

    return decorator
