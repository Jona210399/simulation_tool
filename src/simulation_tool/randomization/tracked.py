import json
from functools import wraps
from pathlib import Path
from typing import Callable, TypedDict


class RandomizationEntry(TypedDict):
    method: str
    args: tuple


type VariableEntry = str

RANDOMIZATION_LOG: dict[VariableEntry, RandomizationEntry] = {}


def record(
    variable_name: str,
    method: str,
    args: tuple,
):
    RANDOMIZATION_LOG[variable_name] = {
        "method": method,
        "args": args,
    }


def reset():
    RANDOMIZATION_LOG.clear()


def export_json(path: Path):
    if path.exists():
        return

    with open(path, "w") as f:
        json.dump(RANDOMIZATION_LOG, f, indent=2)


def tracked(method_name: str):
    def decorator(func: Callable):
        @wraps(func)
        def wrapper(variable_name: str, *args, **kwargs):
            value = func(*args, **kwargs)
            record(variable_name, method_name, args)
            return value

        return wrapper

    return decorator
