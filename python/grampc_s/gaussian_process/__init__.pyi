from __future__ import annotations
import numpy

from typing import Annotated, overload
from numpy.typing import NDArray, ArrayLike

from .kernel import StationaryKernel

__all__: list[str] = ['GaussianProcess', 'GaussianProcessData', 'kernel']
class GaussianProcess:
    @overload
    def __init__(self, data: GaussianProcessData, kernel: StationaryKernel, stateDependency: list[bool], controlDependency: list[bool]) -> None:
        ...
    @overload
    def __init__(self, fileName: str, kernel: StationaryKernel, stateDependency: list[bool], controlDependency: list[bool]) -> None:
        ...
    def dmean_du_vec(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"], vec: float) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    def dmean_dx_vec(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"], vec: float) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    def dvar_du_vec(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"], vec: float) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    def dvar_dx_vec(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"], vec: float) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    def mean(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"]) -> float:
        ...
    def variance(self, evaluationPointState: Annotated[NDArray[numpy.float64], "[m, 1]"], evaluationPointControl: Annotated[NDArray[numpy.float64], "[m, 1]"]) -> float:
        ...

class GaussianProcessData:
    @property
    def X(self) -> Annotated[NDArray[numpy.float64], "[m, n]"]:
        ...
    @X.setter
    def X(self, arg0: Annotated[ArrayLike, numpy.float64, "[m, n]"]) -> None:
        ...
    @property
    def Y(self) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    @Y.setter
    def Y(self, arg0: Annotated[ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def dim(self) -> int:
        ...
    @dim.setter
    def dim(self, arg0: int) -> None:
        ...
    @property
    def output_variance(self) -> float:
        ...
    @output_variance.setter
    def output_variance(self, arg0: float) -> None:
        ...
    @property
    def size(self) -> int:
        ...
    @size.setter
    def size(self, arg0: int) -> None:
        ...
