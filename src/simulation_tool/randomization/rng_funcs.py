from typing import TypedDict

import numpy as np
from numpy.random import Generator, SeedSequence

from simulation_tool.randomization.tracked import tracked


class RangeDict(TypedDict):
    low: float
    high: float


def loguniform(*, rng: Generator, low: float = 0.0, high: float = 1.0, size=None):
    return np.exp(rng.uniform(np.log(low), np.log(high), size))


def uniform(*, rng: Generator, low: float = 0.0, high: float = 1.0):
    return rng.uniform(low, high)


def randn(*, rng: Generator):
    return rng.standard_normal()


def randint(*, rng: Generator, low: int = 0, high: int = 1):
    return rng.integers(low, high)


class RNGFunctions:
    def __init__(self, seed: int | SeedSequence):
        self.rng = np.random.default_rng(seed)

    @tracked(is_method=True)
    def loguniform(self, *, low: float = 0.0, high: float = 1.0, size=None):
        return loguniform(rng=self.rng, low=low, high=high, size=size)

    @tracked(is_method=True)
    def uniform(self, *, low: float = 0.0, high: float = 1.0):
        return uniform(rng=self.rng, low=low, high=high)

    @tracked(is_method=True)
    def randn(self):
        return randn(rng=self.rng)

    @tracked(is_method=True)
    def randint(self, *, low: int = 0, high: int = 1):
        return randint(rng=self.rng, low=low, high=high)
