Introduction
============

This manual describes the structure and the application of GRAMPC-S, a framework for stochastic model predictive control of nonlinear continuous-time systems subject to nonlinear probabilistic constraints.
It is based on the toolbox GRAMPC :footcite:p:`Englert.2019`, which solves an optimization problem using the augmented Lagrangian approach and a tailored gradient method.
GRAMPC-S is implemented as C++ code which allows implementation on embedded hardware and can be used in Simulink and dSPACE.

The documentation is outlined as follows.
:ref:`chap:install` describes the installation of **GRAMPC-S** for Windows and Linux. 
:ref:`chap:system` describes the class of optimal control problems (OCP) that are considered.
:ref:`sec:implementation_problem` shows the implementation of an OCP as C++ code.
Several methods are implemented in GRAMPC-S in order to solve the OCP and propagate uncertainties, whose application is described in :ref:`chap:deterministic_approx`.
For a more detailed explanation of these methods, please refer to :footcite:p:`Landgraf.2023` and :footcite:p:`GRAMPC-S`.
:ref:`chap:simulation` describes the realization of closed-loop simulations and the plotting of results in Matlab.
:ref:`chap:example` concludes with an example of the application and parameterization of GRAMPC-S.

.. footbibliography::
    