from __future__ import annotations
import numpy

from typing import Annotated, overload
from numpy.typing import NDArray

from pygrampc._core import ProblemDescription
from .constraints import ChanceConstraintApproximation
from .distributions import Distribution
from .gaussian_process import GaussianProcess
from .transformations import PointTransformation

__all__: list[str] = ['MonteCarloProblemDescription', 'ResamplingGPProblemDescription', 'ResamplingProblemDescription', 'SigmaPointProblemDescription', 'TaylorBaseProblemDescription', 'TaylorProblemDescription']

class MonteCarloProblemDescription(ProblemDescription):
    def __init__(self, problemDescription: ProblemDescription, pointTransformation: PointTransformation) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution, param: Distribution) -> None:
        ...
    @property
    def p0(self) -> float:
        ...
    @property
    def x0(self) -> float:
        ...

class ResamplingGPProblemDescription(ProblemDescription):
    @overload
    def __init__(self, problemDescription: ProblemDescription, constraintApproximation: ChanceConstraintApproximation, pointTransformation: PointTransformation, gaussianProcessVec: list[GaussianProcess], dynamicsIndicesWithGP: list[int]) -> None:
        ...
    @overload
    def __init__(self, problemDescription: ProblemDescription, pointTransformation: PointTransformation, gaussianProcessVec: list[GaussianProcess], dynamicsIndicesWithGP: list[int]) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution, param: Distribution) -> None:
        ...
    @property
    def p0(self) -> float:
        ...
    @property
    def x0(self) -> float:
        ...

class ResamplingProblemDescription(ProblemDescription):
    @overload
    def __init__(self, problemDescription: ProblemDescription, constraintApproximation: ChanceConstraintApproximation, pointTransformation: PointTransformation) -> None:
        ...
    @overload
    def __init__(self, problemDescription: ProblemDescription, constraintApproximation: ChanceConstraintApproximation, pointTransformation: PointTransformation, diffMarixWienerProcess: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @overload
    def __init__(self, problemDescription: ProblemDescription, pointTransformation: PointTransformation) -> None:
        ...
    @overload
    def __init__(self, problemDescription: ProblemDescription, pointTransformation: PointTransformation, diffMarixWienerProcess: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution, param: Distribution) -> None:
        ...
    @property
    def p0(self) -> float:
        ...
    @property
    def x0(self) -> float:
        ...

class SigmaPointProblemDescription(ProblemDescription):
    @overload
    def __init__(self, problemDescription: ProblemDescription, constraintApproximation: ChanceConstraintApproximation, pointTransformation: PointTransformation) -> None:
        ...
    @overload
    def __init__(self, problemDescription: ProblemDescription, pointTransformation: PointTransformation) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution, param: Distribution) -> None:
        ...
    @property
    def p0(self) -> float:
        ...
    @property
    def x0(self) -> float:
        ...

class TaylorBaseProblemDescription(ProblemDescription):
    pass

class TaylorProblemDescription(ProblemDescription):
    @overload
    def __init__(self, arg0: TaylorBaseProblemDescription, constraintApproximation: ChanceConstraintApproximation) -> None:
        ...
    @overload
    def __init__(self, arg0: TaylorBaseProblemDescription, constraintApproximation: ChanceConstraintApproximation, diffMarixWienerProcess: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @overload
    def __init__(self, arg0: TaylorBaseProblemDescription) -> None:
        ...
    @overload
    def __init__(self, arg0: TaylorBaseProblemDescription, arg1: Annotated[NDArray[numpy.float64], "[m, n]", "flags.f_contiguous"]) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution) -> None:
        ...
    @overload
    def compute_x0_and_p0(self, state: Distribution, param: Distribution) -> None:
        ...
    @property
    def p0(self) -> float:
        ...
    @property
    def x0(self) -> float:
        ...
