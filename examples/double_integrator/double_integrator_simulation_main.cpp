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
#include "double_integrator_problem_description.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
    // State distribution
    Vector x0(2);
    x0 << 2, 0;
    Matrix P0 = 1e-6 * Matrix::Identity(2, 2);
    DistributionPtr state = Gaussian(x0, P0);

    Vector diffusionValues(2);
    diffusionValues << 0.001, 0.001;
    Matrix WienerProcessDiffusionMatirx = diffusionValues.asDiagonal();

    // Point transformation that can be used for the uncertainty propagation
    // PointTransformationPtr transform = UT(state->dimension(), state->dimension(), 0.1, 2, 1);
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension(), state->dimension(), 1);
    PointTransformationPtr transform = StirlingSecondOrder(state->dimension(), state->dimension(), 1);

    // Prediction horizon, sampling time, and initial time
    typeRNum Thor = 4;
    typeRNum dt_MPC = 0.01;
    typeRNum dt_simulation = 0.001;
    typeRNum Tsim = 6;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {0, 0};
    std::vector<typeRNum> udes = {0};

    // Input constraints
    std::vector<typeRNum> umax = {1};
    std::vector<typeRNum> umin = {-1};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pCost;
    std::vector<typeRNum> pCon;

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {0};
    std::vector<typeRNum> constraintsAbsTol = {1e-6};

    // Iterative Cost
    pCost.push_back(1e-3); // Coefficient for (u-u_Des)^2
    pCost.push_back(1.0);  // Coefficient for (x1-x1_Des)^2
    pCost.push_back(1.0); // Coefficient for (x2-x2_Des)^2

    // Terminal Cost 
    pCost.push_back(10.0);  // Coefficient for (x1-x1_Des)^2
    pCost.push_back(10.0);  // Coefficient for (x2-x2_Des)^2
    pCost.push_back(1.0);  // Coefficient for T

    // Constraint parameters
    pCon.push_back(-0.8);

    // Chance constraint probability
    Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.84134;

    // Chance constraint approximation
    // ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // nominal problem description
    ProblemDescriptionPtr integratorProblem = ProblemDescriptionPtr(new DoubleIntegratorProblemDescription(pCost, pCon));

    Vector out(6);
    Vector adj = Vector::Ones(20,1);
    // stochastic problem description
    ResamplingProblemDescriptionPtr problem = ResamplingProblem(integratorProblem, constraintApprox, transform, WienerProcessDiffusionMatirx);
   
    // simulator
    SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, integratorProblem, _1, _2, _3, _4, _5);
    Simulator sim(state->mean(), u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);

    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();

    // set initial states and parameters
    problem->compute_x0_and_p0(state);
    solver->setparam_real_vector("x0", problem->x0());
    solver->setparam_real_vector("p0", problem->p0());

    // resize vectors
    xdes.resize(par->Nx, 0.0);
    constraintsAbsTol.resize(par->Nc, 0.0);

    // set parameters
    solver->setparam_real("Thor", Thor);
    solver->setparam_real("Tmax", 6.0);
    solver->setparam_real("Tmin", 2.0);
    solver->setparam_real("dt", dt_MPC);
    solver->setparam_real("t0", 0.0);
    solver->setparam_real_vector("xdes", &xdes[0]);
    solver->setparam_real_vector("u0", &u0[0]);
    solver->setparam_real_vector("udes", &udes[0]);
    solver->setparam_real_vector("umax", &umax[0]);
    solver->setparam_real_vector("umin", &umin[0]);

    // set options
    solver->setopt_int("MaxGradIter", 2);
    solver->setopt_int("MaxMultIter", 2);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "heun");
    solver->setopt_real("PenaltyMin", 1e3);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
    solver->setopt_string("OptimTime", "on");

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
