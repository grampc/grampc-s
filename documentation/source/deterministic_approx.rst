.. _chap:deterministic_approx:

Approximation of the stochastic OCP
===================================

The :ref:`OCP <eq:ocp_test>` cannot be solved with common solvers, because the calculation of the probabilistic constraints requires the exact knowledge of the probability density functions of the predicted states.
However, the exact propagation of the probability distributions for nonlinear functions can only be calculated by solving a partial differential equation.
Therefore, several approximate methods of uncertainty propagation are implemented in GRAMPC-S.
Their application is explained in Section :ref:`sec:implementation_random_variable`.
Since these only allow an evaluation of a finite number of stochastic moments, the probabilistic constraints must be approximated, which is described in Section :ref:`sec:approx_constraints`.

.. _sec:approx_constraints:

Approximation of probabilistic constraints
------------------------------------------

For the computationally efficient solution of the :ref:`OCP <eq:ocp_test>`, the probabilistic constraints must be converted into deterministic form.
One possibility for this is the Chebyshev inequality, which specifies an upper limit for the probability of a deviation of a random variable from its mean based on the variance of the random variable.
The class ``ChebyshevConstraintApproximation`` implements the approximation of a constraint with the Chebyshev inequality.
An object of this class can be created using the constructor

.. code-block:: C++

	ChebyshevConstraintApproximation(const Vector& probabilities)

where ``probabilities`` describes the vector of probabilities :math:`\vm \alpha` in the probabilistic constraints.
If, in addition to the constraints, terminal constraints are to be taken into account, the concatenation of the vectors :math:`\vm \alpha` and :math:`\vm \alpha_T` is passed to the constructor.
The Chebyshev inequality is valid for any probability density distribution with finite variance, so it is correspondingly conservative.
A less conservative approximation of the probabilistic constraints can be realized if additional information about the distributions of the constraints is available.
The approximation ``SymmetricConstraintApproximation`` assumes that the probability density distribution of the constraints is symmetric.
The approximation can be even less conservative if it is assumed that the constraint distribution is Gaussian.
This type of approximation is implemented in GRAMPC-S in the class ``GaussianConstraintApproximation``.
Objects of these classes can be created with

.. code-block:: C++

	SymmetricConstraintApproximation(const Vector& probabilities)
	GaussianConstraintApproximation(const Vector& probabilities)

.. _sec:stochastic_probabiliy_description:

Stochastic problem descriptions and uncertainty propagation
-----------------------------------------------------------

The system dynamics :math:`\vm f` are given in the form of a stochastic nonlinear differential equation.
Common MPC solvers such as GRAMPC are unable to predict the probability distributions of the states for this class of differential equations.
GRAMPC-S therefore implements several approximate methods of uncertainty propagation that transform the stochastic differential equation into an ordinary differential equation.
These can be divided into two categories: Approximations of the nonlinear function :math:`\vm f` and approximations of the stochastic moments by point-based representations.
The linearization of the system function is explained in Section :ref:`subsec:taylor`, the propagation of uncertainties based on sampling points is explained in Section :ref:`subsec:sigma_point`.

.. _subsec:taylor:

First-order Taylor series approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The consideration of nonlinear functions :math:`\vm f` complicates the calculation of future state distributions. 
This can be easily avoided by linearizing the system dynamics :math:`\vm f` around the current state and calculating the uncertainty propagation with the linearized dynamics.
GRAMPC-S provides the class ``TaylorProblemDescription`` for this purpose.
The constructor

.. code-block:: C++

	TaylorProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef diffMatrixWienerProcess)

receives a pointer ``problemDescription`` to the problem description that implements the system dynamics, the cost function and the constraints.
Since the linearization of the system dynamics is based on the gradients of the system dynamics :math:`\vm f` and GRAMPC in turn solves the optimization problem based on gradients, second-order derivatives are required in addition to the gradients of the functions, as described in Section :ref:`sec:implementation_problem`.
Optionally, the type of constraint approximation and the diffusion matrix of the Wiener process :math:`\vm \Sigma` are passed by the arguments ``constraint-Approximation`` and ``diffMatrixWienerProcess``, respectively.
Gaussian processes are not yet implemented in the class ``TaylorProblemDescription``. 
 
.. _subsec:sigma_point:

Point-based uncertainty propagation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A further way to propagate uncertainties is to approximate probability density distributions by sampling.
The nonlinear function :math:`\vm f` can be evaluated for a set of sampling points and the stochastic moments of the propagated state distribution can be calculated from the results.
This requires two things in GRAMPC-S: A transformation between stochastic moments and sampling points and a stochastic problem description that utilizes this transformation to approximately solve the :ref:`OCP <eq:ocp_test>`. 

The following transformations between stochastic moments and points are implemented in GRAMPC-S, see :footcite:p:`Landgraf.2023` for details:

- Unscented transformation: This is implemented in the class ``UnscentedTransformation`` and requires the scaling parameters :math:`\alpha`, :math:`\beta`, and :math:`\kappa`.
- First-order Stirling's interpolation: The class ``StirlingInterpolationFirstOrder`` implements the first-order Stirling's interpolation for an adjustable step size.
- Second-order Stirling's interpolation: The Stirling's interpolation of second order can be realized with the class ``StirlingInterpolationSecondOrder``, where the step size can be specified as well.
- Gaussian quadrature: The class ``ComposedQuadrature`` generates quadrature points and associated weights for the multivariate Gaussian quadrature by composing it from the quadrature points of one-dimensional Gaussian quadratures. A family of orthogonal polynomials must be known for the respective one-dimensional distribution and the corresponding one-dimensional quadrature rule must be implemented. GRAMPC-S has implementations of the Gauss-Hermite quadrature for Gaussian distributions and the Gauss-Legendre quadrature for uniform distributions. However, it is possible to extend GRAMPC-S with your own quadrature rules by writing a class that inherits from the interface ``QuadratureRule``. In addition to the vector of families of orthogonal polynomials, the constructor receives a vector with the orders of the quadrature rules. If a scalar is passed instead of a vector, this order is used for all quadrature rules.
- Polynomial chaos expansion (PCE): The class ``PCE_Transformation`` contains an implementation of polynomial chaos expansion, which solves the corresponding integrals by Gaussian quadrature and is therefore considered to be a point transformation as well. Its constructor takes the same arguments as the class ``ComposedQuadrature`` and also requires the specification of the maximum polynomial order.
- Monte Carlo sampling: A stochastic sampling of a probability density distribution is implemented in the class ``MonteCarloTransformation``. The number of sampling points must be defined and a random number generator must be passed to the constructor.

It is possible to implement your own point transformations and integrate them into **GRAMPC-S**.
These must be derived classes of the abstract class ``PointTransformation``.
It should be noted that in addition to the calculation of the stochastic moments and points, the associated gradients must also be implemented, because the OCP is solved with a gradient-based approach.
Some of the above point transformations have a second constructor that receives an argument ``considerUncertain``.
This is a vector of boolean values that sets for each variable whether it should be assumed to be uncertain or not. If the respective value is set to ``false``, it is considered as a deterministic variable with the value of the mean of the distribution.
Uncertainties whose impact is negligible can be ignored in this way, which reduces the number of points and thus the required computation time of the controller.
The same can be achieved for Gaussian quadrature and PCE by setting the corresponding entry in vector ``quadratureOrder`` to 1.

The above transformations can be used to represent probability density distributions by points and to calculate stochastic moments from transformed points.
These transformations are used in several problem descriptions that convert the stochastic OCP into a deterministic OCP.
The problem descriptions differ in particular in the internal representation of the states and in whether a Wiener process or a Gaussian process can be taken into account.
An available implementation is the ``SigmaPointProblemDescription``, which initially represents the random variables by points and predicts these points up to the prediction horizon.
The implementation of the system dynamics, the cost function and the constraints are passed to the constructor

.. code-block:: C++

	SigmaPointProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, PointTransformationPtr pointTransformation)

as a pointer to a problem description.
In addition, the type of constraint approximation and the transformation between points and stochastic moments are passed.
The ``SigmaPointProblemDescription`` is computationally efficient, as the computation of points to represent the random variables is only performed once per time step.
However, it is not possible to consider Gaussian processes or an additive Wiener process.
To take these into account, it is not sufficient to initially convert the random variables into sampling points.
Instead, differential equations for the mean and the covariance matrix of the states must be derived, which is implemented in the classes ``ResamplingProblemDescription`` and ``ResamplingGPProblemDescription``.
The integration of these differential equations requires a transformation of stochastic moments into sampling points at each discretization point of the numerical integration.
The system dynamics :math:`\vm f` is evaluated at these points and the time derivative of the mean and covariance matrix of the states is calculated from the results.
An additive Wiener process can be taken into account by the class ``ResamplingProblemDescription`` by passing the diffusion matrix of the Wiener process to the constructor

.. code-block:: C++

    ResamplingProblemDescription(..., MatrixConstRef diffMatrixWienerProcess)

A disadvantage of ``ResamplingProblemDescription`` compared to ``SigmaPointProblemDescription`` is the increased computation time, because the sampling points must be calculated more frequently.
Furthermore, the class ``ResamplingGPProblemDescription`` can be used to consider a Gaussian process approximation of an unknown part of the system dynamics.
For this purpose, the constructor 

.. code-block:: C++

    ResamplingGPProblemDescription(..., const std::vector<GaussianProcessPtr>&  gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP)

receives a vector of Gaussian processes ``gaussianProcessVec``, which are created as described in Section :ref:`sec:implementation_random_variable`.
The output values of the Gaussian processes are each one-dimensional and are combined to obtain the vector :math:`\vm g`.
However, as a Gaussian process should not always be considered for all state differential equations, the argument ``dynamicsIndicesWithGP`` can be used to specify which state differential equations the Gaussian processes relate to.
The size of vector ``gaussianProcessVec`` must be the same as the size of vector ``dynamicsIndicesWithGP``.

A further approximation of the stochastic :ref:`OCP <eq:ocp_test>` is implemented in the class ``MonteCarlo-ProblemDescription``, which provides the constructor

.. code-block:: C++

	MonteCarloProblemDescription(ProblemDescriptionPtr problemDescription,  PointTransformationPtr pointTransformation)

As in ``SigmaPointProblemDescription``, the random variables are initially converted into sampling points, but the probabilistic constraints are not approximated.
Instead, all sampling points must fulfill the deterministic constraints.
In combination with Monte Carlo sampling, the confidence of the fulfillment of the probabilistic constraints can be stated.
The confidence increases with an increasing number of sampling points, however, the corresponding computation time increases as well.

.. footbibliography::