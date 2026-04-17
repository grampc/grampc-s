from __future__ import annotations
from collections.abc import Sequence
import enum
import numpy

from numpy.typing import NDArray, ArrayLike
from typing import ClassVar, Annotated, overload

__all__: list[str] = ['HermitePolynomialGenerator', 'LegendrePolynomialGenerator', 'MultivariatePolynomial', 'OrthogonalPolynomialGenerator', 'Polynomial', 'PolynomialFamily']

class PolynomialFamily(enum.Enum):
    HERMITE: ClassVar[PolynomialFamily]  # value = <PolynomialFamily.HERMITE: 1>
    LEGENDRE: ClassVar[PolynomialFamily]  # value = <PolynomialFamily.LEGENDRE: 2>
    NONE: ClassVar[PolynomialFamily]  # value = <PolynomialFamily.NONE: 0>
    
class Polynomial:
    @overload
    def __init__(self) -> None:
        ...
    @overload
    def __init__(self, coefficients: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    def add_polynomial(self, polynomial: Polynomial) -> None:
        ...
    def evaluate(self, evaluationPoint: float) -> float:
        ...
    def get_coefficient(self, index: int) -> float:
        ...
    def gradient(self, evaluationPoint: float) -> float:
        ...
    def hessian(self, evaluationPoint: float) -> float:
        ...
    def multiply_polynomial(self, polynomial: Polynomial) -> None:
        ...
    def multiply_scalar(self, factor: float) -> None:
        ...
    def subtract_polynomial(self, polynomial: Polynomial) -> None:
        ...
    @property
    def coefficients(self) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    @coefficients.setter
    def coefficients(self, coefficients: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def num_coefficients(self) -> int:
        ...

class OrthogonalPolynomialGenerator:
    def get_maximum_order(self) -> int:
        ...
    def get_polynomial(self, order: int) -> Polynomial:
        ...
    def get_squared_norm(self, order: int) -> float:
        ...

class HermitePolynomialGenerator(OrthogonalPolynomialGenerator):
    def __init__(self, maxOrder: int) -> None:
        ...

class LegendrePolynomialGenerator(OrthogonalPolynomialGenerator):
    def __init__(self, maxOrder: int) -> None:
        ...

class MultivariatePolynomial:
    def __init__(self, univariatePolynomials: Sequence[Polynomial], univariateSquaredNorm: Sequence[float]) -> None:
        ...
    def evaluate(self, evaluationPoint: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> float:
        ...
    @property
    def num_variables(self) -> int:
        ...
    @property
    def polynomials(self, univariatePolynomials: Sequence[Polynomial], univariateSquaredNorm: Sequence[float]) -> None:
        ...
    @property
    def squared_norm(self) -> float:
        ...
        
        
