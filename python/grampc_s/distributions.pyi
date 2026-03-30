from __future__ import annotations
import numpy

from . import RandomNumberGenerator
from .polynomials import PolynomialFamily

from numpy.typing import NDArray
from typing import overload, Annotated, SupportsInt, SupportsFloat

__all__: list[str] = ['Distribution', 'GaussianDistribution', 'MultivariateDistribution', 'UnivariateBetaDistribution', 'UnivariateChiSquaredDistribution', 'UnivariateExponentialDistribution', 'UnivariateExtremeValueDistribution', 'UnivariateFDistribution', 'UnivariateGammaDistribution', 'UnivariateLogNormalDistribution', 'UnivariateStudentTDistribution', 'UnivariateUniformDistribution', 'UnivariateWeibullDistribution']

class Distribution:
    @overload
    def __init__(self, dim: SupportsInt) -> None:
        ...
    @overload
    def __init__(self, mean: Annotated[NDArray[numpy.float64], "[m, 1]"], covariance: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"], poly: list[PolynomialFamily]) -> None:
        ...
    @overload
    def __init__(self, mean: Annotated[NDArray[numpy.float64], "[m, 1]"], covariance: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"], covChol: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"], poly: list[PolynomialFamily]) -> None:
        ...
    def sample(self, rng: RandomNumberGenerator) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    def set_mean_and_covariance(self, mean: Annotated[NDArray[numpy.float64], "[m, 1]"], covariance: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @property
    def cov_cholesky(self) -> Annotated[NDArray[numpy.float64], "[m, n]"]:
        ...
    @property
    def covariance(self) -> Annotated[NDArray[numpy.float64], "[m, n]"]:
        ...
    @property
    def dimension(self) -> int:
        ...
    @property
    def mean(self) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def polynomial_families(self) -> list[PolynomialFamily]:
        ...

class GaussianDistribution(Distribution):
    @overload
    def __init__(self, mean: Annotated[NDArray[numpy.float64], "[m, 1]"], covariance: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @overload
    def __init__(self, mean: SupportsFloat, variance: SupportsFloat) -> None:
        ...

class MultivariateDistribution(Distribution):
    def __init__(self, distributions: list[Distribution]) -> None:
        ...
    def get_distribution(self, index: SupportsInt) -> Distribution:
        ...
    def replace_distribution(self, index: SupportsInt, newDist: Distribution) -> None:
        ...

class UnivariateBetaDistribution(Distribution):
    def __init__(self, p: SupportsFloat, q: SupportsFloat) -> None:
        ...
        
class UnivariateChiSquaredDistribution(Distribution):
    def __init__(self, n: SupportsFloat) -> None:
        ...
        
class UnivariateExponentialDistribution(Distribution):
    def __init__(self, Lambda: SupportsFloat) -> None:
        ...
        
class UnivariateExtremeValueDistribution(Distribution):
    def __init__(self, a: SupportsFloat, b: SupportsFloat) -> None:
        ...
        
class UnivariateFDistribution(Distribution):
    def __init__(self, m: SupportsFloat, n: SupportsFloat) -> None:
        ...
        
class UnivariateGammaDistribution(Distribution):
    def __init__(self, alpha: SupportsFloat, gamma: SupportsFloat) -> None:
        ...
        
class UnivariateLogNormalDistribution(Distribution):
    def __init__(self, mu: SupportsFloat, sigma: SupportsFloat) -> None:
        ...

class UnivariatePieceWiseConstantDistribution(Distribution):
    def __init__(self, intervalLimits: list[SupportsFloat], intervalProbabilityDensity: list[SupportsFloat]) -> None:
        ...
        
class UnivariateStudentTDistribution(Distribution):
    def __init__(self, nu: SupportsFloat, mu: SupportsFloat = 0.0, sigma: SupportsFloat = 1.0) -> None:
        ...
        
class UnivariateUniformDistribution(Distribution):
    def __init__(self, lowerBound: SupportsFloat, upperBound: SupportsFloat) -> None:
        ...
        
class UnivariateWeibullDistribution(Distribution):
    def __init__(self, a: SupportsFloat, b: SupportsFloat) -> None:
        ...
