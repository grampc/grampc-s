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
#include "crane_2d_problem_description.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
    // State  and parameter distribution
    Vector x0(6);
    x0 << -2.0, 0.0, 2.0, 0.0, 0.0, 0.0;
    Matrix P0 = 1e-4 * Matrix::Identity(6, 6);
    DistributionPtr state = Gaussian(x0, P0);

    // Parameter for the propagation methods
    typeInt numberOfPointsMC = 100;
    RandomNumberGenerator rng = RandomNumberGenerator(0);

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = MonteCarlo(state->dimension(), state->dimension(), numberOfPointsMC, rng);

    // Prediction horizon, sampling time
    typeRNum Thor = 2;
    typeRNum dt_MPC = 1e-2;
    typeRNum dt_simulation = 1e-3;
    typeRNum Tsim = 12.5;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {2.0, 0.0, 2.0, 0.0, 0.0, 0.0};
    std::vector<typeRNum> udes = {0.0, 0.0};

    // Input constraints
    std::vector<typeRNum> umax = {2.0,2.0};
    std::vector<typeRNum> umin = {-2.0, -2.0};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pCost = {1.0, 2.0, 2.0, 1.0, 1.0, 4.0, 0.05, 0.05};
    std::vector<typeRNum> pCon = {0.2, 1.25, 0.3};

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {0.0, 0.0};
    std::vector<typeRNum> constraintsAbsTol = {1e-4, 1e-3, 1e-3};

    // Chance constraint probability
    Vector vec_chance_constraint(3);
    vec_chance_constraint << 0.95, 0.95, 0.95;

    // Chance constraint approximation
    // ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    ProblemDescriptionPtr craneProblem = ProblemDescriptionPtr(new Crane2DProblemDescription(pCost, pCon));

    // configure stochastic problem description
    MonteCarloProblemDescriptionPtr problem = MonteCarloProblem(craneProblem, transform);

    // simulator
    SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, craneProblem, _1, _2, _3, _4, _5);
    Simulator sim(state->mean(), u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);
    
    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();

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
    solver->setopt_int("MaxGradIter", 3);
    solver->setopt_int("MaxMultIter", 1);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "heun");
    solver->setopt_real("PenaltyMin", 10);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
    solver->setopt_string("TerminalCost", "off");

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