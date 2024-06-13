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
#include "vehicle_problem_description.hpp"

using namespace grampc;

int main()
{
    // State  and parameter distribution
    Vector x0(5);
    x0 << 10.0, 1.5, 0.0, 0.0, 17.0;
    Vector measurementNoiseVariance(5);
    measurementNoiseVariance << 0.01, 0.01, 1e-5, 1e-5, 1e-3;
    Matrix P0 = measurementNoiseVariance.asDiagonal();
    DistributionPtr state = Gaussian(x0, P0);

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = UT(state->dimension(), state->dimension(), 1, 2, 1);
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = StirlingSecondOrder(state->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = Quadrature(state->dimension(), state->dimension(), state->polynomialFamily(), 3);

    // Prediction horizon, sampling time
    typeRNum Thor = 1;
    typeRNum dt = 0.01;

    std::vector<typeRNum> xdes = {0.0, 1.5, 0.0, 0.0, 0.0};
    std::vector<typeRNum> udes = {0, 0};

    // Input constraints
    std::vector<typeRNum> umax = {0.5, 5.0};
    std::vector<typeRNum> umin = {-0.5,-5.0};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pSys = {3, 50};
    std::vector<typeRNum> pCost = {0.0, 1.0, 1.0, 1.0, 100.0, 0.0, 1.0, 1.0, 1.0, 100.0, 0.01, 0.01};

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {0.0, 0.0};
    std::vector<typeRNum> constraintsAbsTol = {1e-2};

    // Chance constraint probability
    Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.95;

    // Chance constraint approximation
    // ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    ProblemDescriptionPtr vehicleProblem = ProblemDescriptionPtr(new VehicleProblemDescription(pSys, pCost));

    // configure stochastic problem description
    SigmaPointProblemDescriptionPtr problem = SigmaPointProblem(vehicleProblem, constraintApprox, transform);

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
    solver->setparam_real("dt", dt);
    solver->setparam_real_vector("xdes", &xdes[0]);
    solver->setparam_real_vector("u0", &u0[0]);
    solver->setparam_real_vector("udes", &udes[0]);
    solver->setparam_real_vector("umax", &umax[0]);
    solver->setparam_real_vector("umin", &umin[0]);

    // set options
    solver->setopt_int("MaxGradIter", 5);
    solver->setopt_int("MaxMultIter", 1);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "heun");
    solver->setopt_real("PenaltyMin", 1e2);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);

    // computation time
    typeRNum compTime_ms = 0.0;

    // run gradient iterations and save computation time
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    solver->run();
    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    compTime_ms += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() / 1000.0;

    // write simulation data to files
    writeTrajectoriesToFile(solver, state->dimension());

    // Print computation time
    std::cout << "Computation time: " << compTime_ms << " ms" << std::endl;
    std::cout << "Program finished." << std::endl;
    return 0;
}