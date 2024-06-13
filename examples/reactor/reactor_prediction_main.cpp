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
#include "reactor_problem_description.hpp"

using namespace grampc;

int main()
{
    // Parameter for the propagation methods
    typeInt numberOfPoints = 100;
    RandomNumberGenerator rng = RandomNumberGenerator(0);

    // Initial state mean and covariacne 
    Vector x0(2);
    Matrix P0(2,2);
    x0 << 0.8, 0.01;
    P0 << 1e-5, 0, 0, 1e-5;
    DistributionPtr state = Gaussian(x0, P0);

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = PCE(state->dimension(), state->dimension(), state->polynomialFamily(), 2, 3);

    // Prediction horizon, sampling time, and initial time
    typeRNum Thor = 0.01;
    typeRNum dt = 1e-4;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {0.215, 0.09};
    std::vector<typeRNum> udes = {19.6};

    // Input constraints
    std::vector<typeRNum> umax = {100};
    std::vector<typeRNum> umin = {10};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pSys;
    std::vector<typeRNum> pCost;
    std::vector<typeRNum> pCon;

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {19.6};
    std::vector<typeRNum> constraintsAbsTol = {1e-6};

    // System parameters
    pSys.push_back(50);
    pSys.push_back(100);
    pSys.push_back(100);

    // Iterative Cost Q
    pCost.push_back(1.0);  // Coefficient for (x1-x1_Des)^2
    pCost.push_back(1.0); // Coefficient for (x2-x2_Des)^2

    // End cost P
    pCost.push_back(0.0); // Coefficient for (x1-x1_Des)^2
    pCost.push_back(0.0); // Coefficient for (x2-x2_Des)^2

    // Input cost R
    pCost.push_back(2); // Coefficient for (u-u_Des)^2

    // Constraint parameters
    pCon.push_back(0.14);

    // Chance constraint probability
    Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.95;

    // Chance constraint approximation
    ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    ProblemDescriptionPtr reactor_problem = ProblemDescriptionPtr(new ReactorProblemDescription(pSys, pCost, pCon));

    // configure stochastic problem description
    SigmaPointProblemDescriptionPtr problem = SigmaPointProblem(reactor_problem, constraintApprox, transform);

    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();

    // set initial states and parameters depending on the propagation method
    problem->compute_x0_and_p0(state);
    solver->setparam_real_vector("x0", problem->x0());

    // resize vector
    xdes.resize(par->Nx, 0.0);
    constraintsAbsTol.resize(par->Nc, 0.0);

    // set parameters
    solver->setparam_real("Thor", Thor);
    solver->setparam_real("dt", dt);
    solver->setparam_real_vector("xdes", &xdes[0]);
    solver->setparam_real_vector("u0", &u0[0]);
    solver->setparam_real_vector("udes", &udes[0]);
    solver->setparam_real_vector("umax", &umax[0]);
    solver->setparam_real_vector("umin", &umin[0]);

    // set options
    solver->setopt_int("MaxGradIter", 3);
    solver->setopt_int("MaxMultIter", 3);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "heun");
    solver->setopt_real("PenaltyMin", 6e8);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
    solver->setopt_string("InequalityConstraints", "on");
    solver->setopt_string("EqualityConstraints", "off");

    // open output files
    std::ofstream tout("tvec.txt");
    std::ofstream xout("xvec.txt");
    std::ofstream adjout("adjvec.txt");
    std::ofstream uout("uvec.txt");
    std::ofstream dimOut("dim.txt");
    std::ofstream conOut("constr.txt");

    // computation time
    typeRNum compTime_ms = 0.0;

    // number of states, inputs, and constraints
    typeInt Nx = solver->getParameters()->Nx;
    typeInt Nu = solver->getParameters()->Nu;
    typeInt Nc = solver->getParameters()->Nc;

    // run gradient iterations and save computation time
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    solver->run();
    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    compTime_ms += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() / 1000.0;

    writeTrajectoriesToFile(solver, state->dimension());

    // Print computation time
    std::cout << "Computation time: " << compTime_ms << " ms" << std::endl;
    std::cout << "Program finished" << std::endl;
    return 0;
}
