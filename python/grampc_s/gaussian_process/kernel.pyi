from __future__ import annotations
import numpy

from typing import Annotated
from numpy.typing import NDArray

__all__: list[str] = ['KernelProduct', 'KernelSum', 'LocallyPeriodicKernel', 'Matern32Kernel', 'Matern52Kernel', 'PeriodicKernel', 'SquaredExponentialKernel', 'StationaryKernel']

class StationaryKernel:
    def evaluate(self, tau: Annotated[NDArray[numpy.float64], "[m, 1]"]) -> float:
        ...
    def gradient(self, out: Annotated[NDArray[numpy.float64], "[m, 1]", "flags.writeable"], tauEval: Annotated[NDArray[numpy.float64], "[m, 1]"], derivativeIndices: Annotated[NDArray[numpy.int32], "[m, 1]"]) -> None:
        ...
    @property
    def dim(self) -> int:
        ...

class KernelProduct(StationaryKernel):
    def __init__(self, kernel1: StationaryKernel, kernel2: StationaryKernel) -> None:
        ...

class KernelSum(StationaryKernel):
    def __init__(self, kernels: list[StationaryKernel]) -> None:
        ...

class LocallyPeriodicKernel(StationaryKernel):
    def __init__(self, dimInput: int, sigma: float, lengthScale: list[float], period: list[float]) -> None:
        ...

class Matern32Kernel(StationaryKernel):
    def __init__(self, dimInput: int, sigma: float, lengthScale: list[float]) -> None:
        ...

class Matern52Kernel(StationaryKernel):
    def __init__(self, dimInput: int, sigma: float, lengthScale: list[float]) -> None:
        ...

class PeriodicKernel(StationaryKernel):
    def __init__(self, dimInput: int, sigma: float, lengthScale: list[float], period: list[float]) -> None:
        ...

class SquaredExponentialKernel(StationaryKernel):
    def __init__(self, dimInput: int, sigma: float, lengthScale: list[float]) -> None:
        ...

