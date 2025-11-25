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
#include "problem_description_TEMPLATE.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
    // State  and parameter distribution
    Vector x0(1);
    x0 << 0;
    Matrix P0 = 1e-6 * Matrix::Identity(1, 1);
    DistributionPtr state = Gaussian(x0, P0);
    MultiDistributionPtr param = MultiDist({Gaussian(1.0, 1e-3),
                                            Uniform(-1, 1)});

    Vector WienerProcessVariance(1);
    WienerProcessVariance << 0.001;
    Matrix covProcessNoise = WienerProcessVariance.asDiagonal();

    // Parameter for the propagation methods
    typeInt numberOfPointsMC = 100;
    RandomNumberGenerator rng = RandomNumberGenerator(0);

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = UT(state->dimension() + param->dimension(), state->dimension(), 1, 2, 1);
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension() + param->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = StirlingSecondOrder(state->dimension() + param->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = MonteCarlo(state->dimension() + param->dimension(), state->dimension(), numberOfPointsMC, rng);

    // Prediction horizon, sampling time
    typeRNum Thor = 1;
    typeRNum dt_MPC = 0.01;
    typeRNum dt_simulation = 0.001;
    typeRNum Tsim = 3;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {   };
    std::vector<typeRNum> udes = {   };

    // Input constraints
    std::vector<typeRNum> umax = {   };
    std::vector<typeRNum> umin = {   };

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pSys = {   };
    std::vector<typeRNum> pCost = {   };
    std::vector<typeRNum> pCon = {   };

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {   };
    std::vector<typeRNum> constraintsAbsTol = {   };

    // Chance constraint probability
    Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.95;

    // Chance constraint approximation
    // ChanceConstraintApproximationPtr constraintApprox = Chebyshev(vec_chance_constraint);
    // ChanceConstraintApproximationPtr constraintApprox = Symmetric(vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    ProblemDescriptionPtr templateProblem = ProblemDescriptionPtr(new TemplateProblemDescription(pSys, pCost, pCon));

    // configure stochastic problem description
    // ResamplingProblemDescriptionPtr problem = ResamplingProblem(templateProblem, constraintApprox, transform, covProcessNoise);
    // MonteCarloProblemDescriptionPtr problem = MonteCarloProblem(templateProblem, transform);
    SigmaPointProblemDescriptionPtr problem = SigmaPointProblem(templateProblem, constraintApprox, transform);
    
    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();
    
    // simulator
    SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, templateProblem, _1, _2, _3, _4, _5, par);
    Simulator sim(state->mean(), param->mean(), u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);

    // set initial states and parameters depending on the propagation method
    problem->compute_x0_and_p0(state, param);
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
    solver->setopt_int("MaxMultIter", 2);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "erk2");
    solver->setopt_real("PenaltyMin", 1e6);
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
        x0 = sim.simulate(solver->getSolution()->unext);
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