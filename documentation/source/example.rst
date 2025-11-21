.. _chap:example:

Example: Continuous stirred tank reactor
========================================

This section shows the implementation of a control of a continuous stirred tank reactor with GRAMPC-S, which can be found in the folder ``examples/reactor``. A simplified and normalized version of the model from :footcite:p:`KlattEngell:95`, :footcite:p:`ROTHFUSS19961433` is used, which is :footcite:p:`GRAICHEN20041123` 

.. math::
    :label: eq:reactor_dyn

	\dot c_A &= -k_1 c_A - k_3 c_A^2 + (1 - c_A) u

	\dot c_B &= k_1 c_A - k_2 c_B - c_B u \,,

where :math:`c_A` is the normalized concentration of the educt :math:`A`, :math:`c_B` is the normalized concentration of product :math:`B`, :math:`u` is the volume flow rate, and

.. math::

    k_1 = 50\text{h}^-1\,, \qquad k_2 = 100\text{h}^-1\,, \qquad k_3 = 100\text{h}^-1

are parameters of the process.
The optimization problem 

.. math:: 
    :label: eq:ocp_reactor

    \min_{u} \,\,\quad & J\left( u; \, \vm x_0 \right) = \mathbb{E} \left[\int\limits_0^T \!\Delta \vm x \trans \vm Q \Delta \vm x  + 0.1 \Delta u^2 \,\, \text{d}t \right] 

    \operatorname{s.t.} \quad & \dot{\vm{x}} = \vm f(\vm{x}, \, u)\,,\quad \vm{x}(0) = \vm x_0

    & \mathbb{P} \left[ c_B \leq 0.14 \right] \geq 0.95

    & u \in \left[ 10 ,\, 100 \right]

should be solved, where the system dynamics are given by :math:numref:`eq:reactor_dyn` and :math:`\Delta \vm x = \vm x - \vm x_\text{des}` and :math:`\Delta u = u - u_\text{des}` denote the deviation of the states and the control input from the corresponding desired values, with the weight matrix :math:`\vm Q = \mathrm{diag}([10^2, 10^2])`.
It is assumed that the current state is not exactly known due to measurement noise, which makes it necessary to predict the uncertainties.
These are required to evaluate the probabilistic constraint for the concentration :math:`c_B`.
The problem description will now be implemented for this OCP.
For this purpose, a class ``ReactorProblemDescription`` is created, which provides all the necessary functions from Section :ref:`sec:implementation_problem`.
First, the constructor is implemented:

.. code-block:: C++

    ReactorProblemDescription::ReactorProblemDescription(const std::vector<typeRNum>& pSys,
                                                        const std::vector<typeRNum>& pCost, 
                                                        const std::vector<typeRNum>& pCon)
     : pSys_(pSys), pCost_(pCost), pCon_(pCon)
    {}

The constructor receives vectors containing the system parameters, the coefficients for the cost function and for the constraint and saves them in member variables.
Next, the number of states, the number of control inputs, the number of parameters, the number of constraints and the number of terminal constraints of the OCP must be specified.
This is realized in the function ``ocp_dim``:

.. code-block:: C++

    void ReactorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = 2;  *Nu = 1;  *Np = 0;  *Ng = 0;
        *Nh = 1;  *NgT = 0;  *NhT = 0; 
    }

Please note that the generic data type ``typeInt`` for integer values is used here which is defined by GRAMPC in ``grampc_macro.h``.
The generic data types ``typeRNum`` and ``ctypeRNum`` for floating point numbers are also defined there. 
The system dynamics in :math:numref:`eq:reactor_dyn` are implemented in the function ``ffct``:

.. code-block:: C++

    void ReactorProblemDescription::ffct(VectorRef out, ctypeRNum t,  VectorConstRef x, VectorConstRef u, VectorConstRef p)
    {
        out[0] = -pSys_[0] * x[0] - pSys_[2] * x[0] * x[0] + (1 - x[0]) * u[0];
        out[1] = pSys_[0] * x[0] - pSys_[1] * x[1] - x[1] * u[0]; 
    }

Here, the member variable ``pSys_`` is used, which was set in the constructor.
Next, the cost function must be implemented in the function ``lfct``:

.. code-block:: C++

	void ReactorProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x,  VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        out[0] = pCost_[2] * (x[0] - xdes[0]) * (x[0] - xdes[0]) +
        pCost_[3] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +
        pCost_[4] * (u[0] - udes[0]) * (u[0] - udes[0]);
    }

Again, a member variable is accessed, which was set in the constructor.
If terminal constraints need to be taken into account in addition to the iterative costs, these must be implemented in the function ``Vfct``.
The inequality constraint is implemented in the function ``hfct``:

.. code-block:: C++

    void ReactorProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x,  VectorConstRef u, VectorConstRef p)
    { 
        out[0] = x[1] - pCon_[0];  
    }

The underlying solver GRAMPC is gradient-based, so some derivatives with respect to :math:`\vm x` and :math:`u` must be implemented for the above functions.
As described in Section :ref:`sec:implementation_problem`, the matrix product of the Jacobian matrix and an arbitrary vector, which is passed as a pointer, is returned for the derivatives of the system dynamics and the constraint:

.. code-block:: C++

    void ReactorProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t,  VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
    {
        out[0] = (-pSys_[0] - pSys_[2] * 2 * x[0] - u[0]) * adj[0] + pSys_[0] * adj[1];
        out[1] = (-pSys_[1] - u[0]) * adj[1];
    }

    void ReactorProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t,  VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
    {
        out[0] = (1 - x[0]) * adj[0] + (-x[1]) * adj[1];
    }

    void ReactorProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x,  VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        out[0] = 2 * pCost_[2] * (x[0] - xdes[0]);
        out[1] = 2 * pCost_[3] * (x[1] - xdes[1]);
    }

    void ReactorProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x,  VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        out[0] = 2 * pCost_[4] * (u[0] - udes[0]);
    }

    void ReactorProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t,  VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
    {
        out[0] = 0.0;
        out[1] = vec[0];
    }

    void ReactorProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t,  VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
    {
        out[0] = 0.0;
    }

All functions of the OCP are now defined.
The next step is to define the distribution of states, create and parameterize the solver and implement the MPC loop.
This is realized in the function ``main``, which is located in the file ``reactor_simulation_main.cpp`` for this example.
First, the initial state :math:`\vm x_0` is created as a Gaussian distribution with mean :math:`\mathbb{E}[\vm x_0] = [0.7,\, 0.01]\trans` and covariance matrix :math:`P = 10^{-5} \; \vm I`:

.. code-block:: C++

    Vector x0(2); 
    Matrix P0(2,2);
    x0 << 0.7, 0.01;
    P0 << 1e-5, 0, 0, 1e-5;
    DistributionPtr state = Gaussian(x0, P0);

In principle, all methods from Section :ref:`chap:deterministic_approx` can be used to propagate the uncertainties.
In this example we use polynomial chaos expansion with a maximum polynomial order of 2 and solve the occurring integrals using Gaussian quadrature with a grid of 3 points per dimension:

.. code-block:: C++

    PointTransformationPtr transform = PCE(2, 2, state->polynomialFamily(), 2, 3);

Next, the approximation of the probabilistic constraint is determined by the Chebyshev inequality:

.. code-block:: C++

	Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.95;
    ChanceConstraintApproximationPtr constraintApprox =  Chebyshev(vec_chance_constraint);

The class ``ReactorProblemDescription`` has already been implemented in the previous steps.
An object of this class is now created and passed to a ``SigmaPointProblemDescription``:

.. code-block:: C++

    ProblemDescriptionPtr reactorProblem = ProblemDescriptionPtr( new ReactorProblemDescription(pSys, pCost, pCon)); 
    SigmaPointProblemDescriptionPtr problem = SigmaPointProblem( reactorProblem, constraintApprox, transform);

The ``SigmaPointProblemDescription`` represents the approximation of the stochastic OCP, which can be solved by GRAMPC.
To this end, a solver must be created and the pointer to the problem description must be passed to it:

.. code-block:: C++

	GrampcPtr solver = Solver(problem);

Next, the solver must be parameterized, which is skipped here for reasons of space, but can be found in the file ``reactor_simulation_main.cpp``.
Finally, the MPC loop is implemented.
In each time step, the solver is called, the system is simulated and the current state and time are passed to the solver:

.. code-block:: C++

    for(typeRNum t = 0; t < Tsim; t += dt_MPC)
    {
        solver->run();

        //  simulate system
        Eigen::Map<const Vector> uMap(solver->getSolution()->unext, u0.size());
        x0 = sim.simulate(uMap);
        state->setMeanAndCovariance(x0, P0);

        // set current state and time
        problem->compute_x0_and_p0(state);
        solver->setparam_real_vector("x0", problem->x0());
        solver->setparam_real("t0", t);
    }

:numref:`fig:plot_results_reactor` shows the results of the simulation for desired states :math:`\vm x_\text{des} = [0.215, \,0.09]\trans` and a desired control input :math:`u_\text{des} = 19.6`.
First, the concentration :math:`c_A` decreases and the concentration :math:`c_B` increases until the constraint becomes active.
The deviation between the trajectory of :math:`c_B` and the value of 0.14 results from the constraint tightening with the Chebyshev inequality due to the state uncertainty.

.. figure:: figures/reactor_example/plot_reactor_results.*
    :name: fig:plot_results_reactor

    Results of the stochastic MPC simulation for the chemical reactor.

.. footbibliography::
    