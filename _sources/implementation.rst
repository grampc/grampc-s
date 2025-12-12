.. _chap:implementation:

Implementation
====================================

The optimization problem must be implemented as C++ code in order to be used in GRAMPC-S.
All functions that appear in the OCP must be defined in a problem description class, which is described in Section :ref:`sec:implementation_problem`.
In addition, random variables such as initial states and parameters must be defined in each time step and, optionally, the Gaussian process :math:`\vm g` must be specified, which is described in Section :ref:`sec:implementation_random_variable`.


.. _sec:implementation_problem:

Implementation of the problem description
-----------------------------------------

The structure of the :ref:`OCP <eq:ocp_test>` consisting of the system dynamics, costs and constraints must be implemented as methods of a class that inherits from the abstract class ``ProblemDescription``. The **examples** folder contains implementations of several examples as well as the template ``problem_description_TEMPLATE.hpp`` for the header file and the template ``problem_description_TEMPLATE.cpp`` for the source code. The following functions must be implemented:

- ``ocp_dim``: Definition of the dimensions of the :ref:`OCP <eq:ocp_test>`, such as the number of states, the number of control inputs, the number of constraints, etc.
- ``ffct``: The system dynamics :math:`\vm f`. 
- ``dfdx_vec``: The matrix vector product :math:`\left( \frac{\partial \vm f}{\partial \vm x}\right)\trans \vm v` of the Jacobian :math:`\frac{\partial \vm f}{\partial \vm x}` with an arbitrary vector :math:`\vm v`.

Note that by using the matrix product of Jacobian matrix and vector, excessive multiplications with zero should be avoided. Depending on the OCP, further optional parts can be implemented:

- ``dfdu_vec``: The matrix product :math:`\left( \frac{\partial \vm f}{\partial \vm u}\right)\trans \vm v` of the Jacobian :math:`\frac{\partial \vm f}{\partial \vm u}` with an arbitrary vector :math:`\vm v`.
- ``lfct``, ``Vfct``: The integral cost :math:`l` and the terminal cost :math:`V`.
- ``dldx``, ``dldu``: The gradient of the integral cost with respect to :math:`\vm  x` and :math:`\vm u`.
- ``dVdx``, ``dVdT``: The gradient of the terminal cost with respect to :math:`\vm  x` and :math:`T`.
- ``hfct``, ``hTfct``: The inequality constraints :math:`\vm h` and :math:`\vm h_T`.
- ``dhdx_vec``, ``dhdu_vec``, ``dhdp_vec``: The matrix vector products :math:`\left( \frac{\partial \vm h}{\partial \vm x}\right)\trans \vm v`, :math:`\left( \frac{\partial \vm h}{\partial \vm u}\right)\trans \vm v`, and :math:`\left( \frac{\partial \vm h}{\partial \vm p}\right)\trans \vm v` of the Jacobians with an arbitrary vector :math:`\vm v`.
- ``dhTdx_vec``, ``dhTdp_vec``, ``dhTdT_vec``: The matrix vector products :math:`\left( \frac{\partial \vm h_T}{\partial \vm x}\right)\trans \vm v`, :math:`\left( \frac{\partial \vm h_T}{\partial \vm T}\right)\trans \vm v`, and :math:`\left( \frac{\partial \vm h_T}{\partial \vm p}\right)\trans \vm v` of the Jacobians with an arbitrary vector :math:`\vm v`.

The above functions are sufficient for most of the propagation methods in GRAMPC-S. However, if you want to use the Taylor series expansion around the current state in order to propagate uncertainties (see :ref:`chap:deterministic_approx` for details), the following functions must also be implemented:

- ``dfdx``, ``dfdp``: Jacobians :math:`\frac{\partial \vm f}{\partial \vm x}` and :math:`\frac{\partial \vm f}{\partial \vm p}`.
- ``dfdxdx``, ``dfdxdp``, ``dfdxdu``, ``dfdpdu``: Hessian matrices :math:`\frac{\partial^2 \vm f}{\partial^2 \vm x}`, :math:`\frac{\partial^2 \vm f}{\partial \vm x \partial \vm p}`, :math:`\frac{\partial^2 \vm f}{\partial \vm x \partial \vm u}`, and :math:`\frac{\partial^2 \vm f}{\partial \vm p \partial \vm u}`.
- ``dhdx``, ``dhdu``: Jacobians :math:`\frac{\partial \vm h}{\partial \vm x}` and :math:`\frac{\partial \vm h}{\partial \vm u}`.
- ``dhdxdx``, ``dhdxdu``: Hessian matrices :math:`\frac{\partial^2 \vm h}{\partial^2 \vm x}` and :math:`\frac{\partial^2 \vm h}{\partial \vm x \partial \vm u}`.
- ``dhTdx``: Jacobian matrix :math:`\frac{\partial \vm h_T}{\partial \vm x}`.
- ``dhTdxdx``, ``dhTdxdT``: Hessian matrices :math:`\frac{\partial^2 \vm h_T}{\partial^2 \vm x}` and :math:`\frac{\partial^2 \vm h_T}{\partial \vm x \partial T}`.

The problem description is usually implemented in a .cpp file. To ensure that this is taken into account in the compilation process, it must be listed in the ``CMakeList.txt`` file. The folder **examples** contains several examples that can be used as references for problem descriptions.

.. _sec:implementation_random_variable:

Implementation of random variables and Gaussian processes
---------------------------------------------------------

The initial states :math:`\vm x_0` are generally random variables, with a given probability density function. :ref:`sec:appendix_pdf` lists all probability density functions that are implemented in GRAMPC-S.
However, it is possible to write your own probability density function as a class, which must be declared as derived from the ``Distribution`` class.
The following example presents the definition of a Gaussian distribution with mean :math:`[1,\, 0]\trans` and covariance matrix :math:`10^{-4} \cdot \vm I`, where :math:`\vm I` is the identity matrix:

.. code-block:: C++

    Vector x0(2); 
    x0 << 1.0, 0.0; 
    Matrix P0 = 1e-4 * Matrix::Identity(2, 2); 
    DistributionPtr state = Gaussian(x0, P0);

The resulting object is a shared pointer to the Gaussian distribution and can be used subsequently. In the closed-loop setting, the distribution must be updated in each time step. The ``Distribution`` class provides the member function ``setMeanAndCovariance`` for this purpose: 

.. code-block:: C++

    Vector x1(2);	
    x1 << 1.5, -1.0; 
    Matrix P1 = 1e-2 * Matrix::Identity(2, 2); 
    state->setMeanAndCovariance(x1, P1);

GRAMPC-S is also able to consider a Gaussian process that takes non-modeled parts of the system dynamics into account.
A Gaussian process is a non-parametric data-based model, see :footcite:p:`bib:rasmussen2005` for details.
First, the kernel and the hyperparameters of the Gaussian process must be defined, selecting a kernel from the list in :ref:`sec:appendix_kernels`.
It is possible to write your own kernel by writing a class that inherits from the ``StationaryKernel`` class.
For example, a squared exponential kernel for a two-dimensional input variable with output variance :math:`\sigma^2 = 1` and length scaling :math:`\vm l = [2,\, 3]\trans` can be defined as follows:

.. code-block:: C++

	StationaryKernelPtr kernel = SquaredExponential(2, 1, {2, 3});

Furthermore, the data points must be passed to the Gaussian process.
There are two options for this: Using the structure ``GaussianProcessData``, which is passed to the Gaussian process, or by reading in a text document in which the input data is stored.
The ``matlab`` folder contains the Matlab function ``writeGaussianProcessData.m``, which can be used to easily create such a text file.
The folder ``examples/double_integrator_GP`` contains an example of a double integrator whose system dynamics are completely represented by two Gaussian processes.
The Matlab script ``double_integrator_data.m`` in this folder shows an example how data can be generated and how the function ``writeGaussianProcessData.m`` should be called.
Once the text document has been generated, a Gaussian process can be created for example with

.. code-block:: C++

	GaussianProcessPtr gp = GP("GP.txt", kernel, {false, true}, {true});

Here, ``GP.txt`` is the name of the text document, ``kernel`` is the pointer to the previously created kernel and the following vectors of type boolean describe the dependency of the Gaussian process on the states and the control input.
In this example, the Gaussian process depends on the second state and the control input.

.. footbibliography::