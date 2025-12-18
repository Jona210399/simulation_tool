import json
from functools import wraps
from pathlib import Path
from typing import Any, Callable, TypedDict


class RandomizationEntry(TypedDict):
    method: str
    args: tuple
    kwargs: dict[str, Any]


type VariableEntry = str

RANDOMIZATION_LOG: dict[VariableEntry, RandomizationEntry] = {}


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


def tracked(preserve_first_arg: bool = False) -> Callable:
    def decorator(func: Callable):
        method_name = func.__qualname__
        if preserve_first_arg:

            @wraps(func)
            def wrapper(self, variable_name: str, *args, **kwargs):
                value = func(self, *args, **kwargs)
                record(variable_name, method_name, args, kwargs)
                return value

            return wrapper

        @wraps(func)
        def wrapper(variable_name: str, *args, **kwargs):
            value = func(*args, **kwargs)
            record(variable_name, method_name, args, kwargs)
            return value

        return wrapper

    return decorator
