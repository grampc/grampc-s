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
    typeRNum dt = 0.005;

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
    solver->setparam_real("dt", dt);
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