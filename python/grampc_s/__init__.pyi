from __future__ import annotations

import numpy

from typing import Annotated, overload
from collections.abc import Callable
from numpy.typing import NDArray

from . import constraints
from . import distributions
from . import gaussian_process
from . import polynomials
from . import problem_descriptions
from . import transformations

__all__: list[str] = ['RandomNumberGenerator', 'Simulator', 'constraints', 'distributions', 'gaussian_process', 'polynomials', 'problem_descriptions', 'transformations']

class RandomNumberGenerator:
    def __init__(self, seed: int) -> None:
        ...

class Simulator:
    @overload
    def __init__(self, initialState: Annotated[NDArray[numpy.float64], "[m, 1]"], param: Annotated[NDArray[numpy.float64], "[m, 1]"], numberOfControls: int, systemFct: Callable[[Annotated[NDArray[numpy.float64], "[m, 1]", "flags.writeable"], float, Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"]], None], integrator: str, diffMatrixWienerProcess: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"], t0: float, dt_MPC: float, dt_simulation: float, writeDateToFile: bool) -> None:
        ...
    @overload
    def __init__(self, initialState: Annotated[NDArray[numpy.float64], "[m, 1]"], numberOfControls: int, systemFct: Callable[[Annotated[NDArray[numpy.float64], "[m, 1]", "flags.writeable"], float, Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"]], None], integrator: str, diffMatrixWienerProcess: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"], t0: float, dt_MPC: float, dt_simulation: float, writeDataToFile: bool) -> None:
        ...
    @overload
    def __init__(self, initialState: Annotated[NDArray[numpy.float64], "[m, 1]"], param: Annotated[NDArray[numpy.float64], "[m, 1]"], numberOfControls: int, systemFct: Callable[[Annotated[NDArray[numpy.float64], "[m, 1]", "flags.writeable"], float, Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"]], None], integrator: str, t0: float, dt_MPC: float, dt_simulation: float, writeDataToFile: bool) -> None:
        ...
    @overload
    def __init__(self, initialState: Annotated[NDArray[numpy.float64], "[m, 1]"], numberOfControls: int, systemFct: Callable[[Annotated[NDArray[numpy.float64], "[m, 1]", "flags.writeable"], float, Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"], Annotated[NDArray[numpy.float64], "[m, 1]"]], None], integrator: str, t0: float, dt_MPC: float, dt_simulation: float, writeDataToFile: bool) -> None:
        ...
    def simulate(self, control: Annotated[NDArray[numpy.float64], "[m, 1]"]) -> Annotated[NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def t(self) -> float:
        ...
