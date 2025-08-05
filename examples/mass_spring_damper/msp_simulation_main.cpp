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
#include "msp_problem_description.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
    // Number of masses
    typeInt numberOfMasses = 5;

    // State  and parameter distribution
    Vector x0 = Vector::Zero(2 * numberOfMasses);
    x0(0) = 1.0;
    Matrix P0 = 1e-9 * Matrix::Identity(2 * numberOfMasses, 2 * numberOfMasses);

    Vector WienerProcessVariance = 1e-2 * Vector::Ones(x0.rows());
    Matrix covProcessNoise = WienerProcessVariance.asDiagonal();

    DistributionPtr state = Gaussian(x0, P0);
    MultiDistributionPtr param = MultiDist({Uniform(0.98, 1.02),
                                            Uniform(0.18, 0.22)});
   
    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = UT(state->dimension() + param->dimension(), state->dimension(), 1, 2, -3);
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension() + param->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = StirlingSecondOrder(state->dimension() + param->dimension(), state->dimension(), 1);

    // Prediction horizon, sampling time
    typeRNum Thor = 2;
    typeRNum dt_MPC = 5e-3;
    typeRNum dt_simulation = 5e-4;
    typeRNum Tsim = 12;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes(2 * numberOfMasses, 0.0);
    std::vector<typeRNum> udes = {0.0, 0.0};

    // Input constraints
    std::vector<typeRNum> umax = {1, 1};
    std::vector<typeRNum> umin = {-1, -1};

    // Declaration of vectors for the problem description
    std::vector<typeRNum> pSys = {1};
    std::vector<typeRNum> pCost(4*numberOfMasses + 2, 1.0);
    pCost[2 * numberOfMasses] = 0.01;
    pCost[2 * numberOfMasses + 1] = 0.01;

    // Initial input and tolerance of the constraints
    std::vector<typeRNum> u0 = {0.0, 0.0};

    // create nominal problem description
    ProblemDescriptionPtr mspProblem = ProblemDescriptionPtr(new MassSpringDamperProblemDescription(numberOfMasses, pSys, pCost));

    // configure stochastic problem description
    SigmaPointProblemDescriptionPtr problem = SigmaPointProblem(mspProblem, transform);
    // TaylorProblemDescriptionPtr problem = TaylorProblem(mspProblem, covProcessNoise);
    // ResamplingProblemDescriptionPtr problem = ResamplingProblem(mspProblem, transform, covProcessNoise);

    // simulator
    SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, mspProblem, _1, _2, _3, _4, _5);
    Simulator sim(state->mean(), param->mean(), u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);
    
    // create solver
    GrampcPtr solver = Solver(problem);
    const typeGRAMPCparam *par = solver->getParameters();

    // set initial states and parameters depending on the propagation method
    problem->compute_x0_and_p0(state, param);
    solver->setparam_real_vector("x0", problem->x0());
    solver->setparam_real_vector("p0", problem->p0());

    // resize vectors
    xdes.resize(par->Nx, 0.0);

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
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("InequalityConstraints", "off");

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