/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


#include "grampc_s.hpp"
#include "inverted_pendulum_problem_description.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
    // definition of pi
    constexpr typeRNum pi = 3.14159265358979323846;

    // State  and parameter distribution
    Vector x0(4);
    Matrix P0 = 1e-6 * Matrix::Identity(4, 4);
    x0 << 0.5, 0, pi, 0.0;

    DistributionPtr state = Gaussian(x0, P0);

    Vector diffusionValues(4);
    diffusionValues << 0, 0.01, 0, 0.01;
    Matrix WienerProcessDiffusionMatirx = diffusionValues.asDiagonal();

    // Prediction horizon, sampling time
    typeRNum Thor = 0.7;
    typeRNum dt_MPC = 0.02;
    typeRNum dt_simulation = 1e-4;
    typeRNum Tsim = 3;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {0, 0, pi, 0};
    std::vector<typeRNum> udes = {0};

    // Input constraints
    std::vector<typeRNum> umax = {5};
    std::vector<typeRNum> umin = {-5};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pSys = {2e-3, 0.2, 0.15, 9.81, 1.3e-3};
    std::vector<typeRNum> pCost = {100, 1, 1, 1, 0, 0, 0, 0, 1e-9};
    std::vector<typeRNum> pCon = {0.8};

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {0};
    std::vector<typeRNum> constraintsAbsTol = {1e-6, 1e-6};

    // Chance constraint probability
    Vector vec_chance_constraint(2);
    vec_chance_constraint << 0.95, 0.95;

    // Chance constraint approximation
    // ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    ProblemDescriptionPtr pendulumProblem = ProblemDescriptionPtr(new InvertedPendulumProblemDescription(pSys, pCost, pCon));

    // configure stochastic problem description
    TaylorProblemDescriptionPtr problem = TaylorProblem(pendulumProblem, constraintApprox, WienerProcessDiffusionMatirx);

    typeInt numberOfPoints = 10000;
    RandomNumberGenerator rng = RandomNumberGenerator(0);

    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();
    
    // simulator
    SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, pendulumProblem, _1, _2, _3, _4, _5, par);
    Simulator sim(state->mean(), (typeInt) u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);

    // set initial states and parameters depending on the propagation method
    problem->compute_x0_and_p0(state);
    solver->setparam_real_vector("x0", problem->x0());
    solver->setparam_real_vector("p0", problem->p0());

    // resize vectors
    xdes.resize(par->Nx, 0.0);
    constraintsAbsTol.resize(par->Nc, 0.0);

    // set parameters
    solver->setparam_real("Thor", Thor);
    solver->setparam_real("dt", dt_MPC);
    solver->setparam_real_vector("xdes", &xdes[0]);
    solver->setparam_real_vector("u0", &u0[0]);
    solver->setparam_real_vector("udes", &udes[0]);
    solver->setparam_real_vector("umax", &umax[0]);
    solver->setparam_real_vector("umin", &umin[0]);

    // set options
    solver->setopt_int("MaxGradIter", 5);
    solver->setopt_int("MaxMultIter", 2);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "erk2");
    solver->setopt_real("LineSearchMax", 10);
    solver->setopt_real("PenaltyMin", 1e3);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);

    // computation time
    typeRNum compTime_ms = 0.0;

    // MPC loop
    for(typeRNum t = 0; t < Tsim; t += dt_MPC)
    {
        // run solver and save computation time
        std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
        solver->run();
        std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
        compTime_ms += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() / 1000.0;

        // simulate system
        Eigen::Map<const Vector> uMap(solver->getSolution()->unext, u0.size());
        x0 = sim.simulate(uMap);
        state->setMeanAndCovariance(x0, P0);

        // set current state and time
        problem->compute_x0_and_p0(state);
        solver->setparam_real_vector("x0", problem->x0());
        solver->setparam_real("t0", t);
    }

    // Print computation time
    std::cout << "Average computation time per sampling step: " << compTime_ms / std::round(Tsim/dt_MPC)<< " ms" << std::endl;
    std::cout << "Program finished." << std::endl;
    return 0;
}