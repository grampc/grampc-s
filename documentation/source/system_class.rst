.. _chap:system:

System class
====================================

GRAMPC-S considers optimal control problems of the form

.. _eq:OCP_test:

.. math::
    :label: eq:OCP

    \min_{\vm u,T} \,\,\quad& J\left( \vm u; \, \vm p, \, \vm x_0 \right) = \mathbb{E} \left[ V(\vm x(T), \, \vm p, \, T) + \int\limits_0^T l(\vm{x}, \vm{u}, \, \vm p, \,t) \, \text{d}t \right]

    \operatorname{s.t.} \quad& \text{d}\vm{x} = \vm f(\vm{x}, \, \vm{u}, \, \vm p, \, t) \, \text{d}t + \vm g(\vm{x}, \, \vm{u}) \, \text{d}t + \vm \Sigma \, \text{d} \vm w \,,\quad \vm{x}(0) = \vm x_0

    & \mathbb{P} \left[h_i(\vm{x}, \, \vm{u}, \ \vm p,\, t) \leq \vm 0 \right] \geq \alpha_i \,,\qquad \quad \hspace{0.48em} i \in \left[ 1,\, n_h \right]

    & \mathbb{P} \left[h_{T,\,j}(\vm{x}(T), \, \vm p, \ T) \leq \vm 0 \right] \geq \alpha_{T,\,j} \,, \quad \, \, j \in \left[ 1,\, n_{h_T} \right]

    & \vm{u} \in \left[ \vm u_{min} ,\, \vm u_{max} \right]

    & T \in \left[T_{min}, \, T_{max}\right]

in the context of model predictive control. The functional of the expected cost :math:`J\left( \vm u; \, \vm p, \, \vm x_0 \right)` consists of the expected value of terminal cost :math:`V(\vm x(T), \, \vm p, \, T)` and integral cost :math:`l(\vm{x}, \vm{u}, \, \vm p, \,t)`, which depend on states :math:`\vm x \in \mathbb{R}^{n_x}`, control inputs :math:`\vm u \in \mathbb{R}^{n_u}`, parameters :math:`\vm p\in \mathbb{R}^{n_p}`, and the prediction horizon :math:`T`.
States and parameters are assumed to be random variables.

The cost functional :math:`J\left( \vm u; \, \vm p, \, \vm x_0 \right)` is minimized subject to the system dynamics, probabilistic constraints, probabilistic terminal constraints, box constraints for :math:`\vm u` and :math:`T`, respectively.
The system dynamics are given in the form of a stochastic differential equation in which :math:`\vm f` describes the part that is known, :math:`\vm g` is a Gaussian process approximation of unknown parts of the system dynamics, and :math:`\vm w` is a multidimensional Wiener process with constant diffusion term :math:`\vm \Sigma`.
The Gaussian process approximation and the Wiener noise are optional components.
If these are not specified, they are not taken into account in order to reduce the computation time.
Based on the system dynamics, the states are predicted starting from the initial state :math:`\vm x_0`.
The vectors :math:`\vm h \in \mathbb{R}^{n_h}` and :math:`\vm h_T \in \mathbb{R}^{n_{h_T}}` denote the inequality constraints and terminal inequality constraints that must be satisfied with probabilities :math:`\vm \alpha \in \mathbb{R}^{n_h}` and :math:`\vm \alpha_T \in \mathbb{R}^{n_{h_T}}`, respectively.  