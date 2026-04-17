from __future__ import annotations
import numpy

from numpy.typing import ArrayLike, NDArray
from typing import Annotated

__all__: list[str] = ['ChanceConstraintApproximation', 'ChebyshevConstraintApproximation', 'GaussianConstraintApproximation', 'SymmetricConstraintApproximation']

class ChanceConstraintApproximation:
    def __init__(self, probabilities: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def probability(self, probabilities: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @probability.setter
    def probability(self) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def z(self) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...

class ChebyshevConstraintApproximation(ChanceConstraintApproximation): ...
class GaussianConstraintApproximation(ChanceConstraintApproximation): ...
class SymmetricConstraintApproximation(ChanceConstraintApproximation): ...