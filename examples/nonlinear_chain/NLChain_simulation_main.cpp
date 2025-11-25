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
#include "NLChain_2_problem_description.hpp"
#include "NLChain_3_problem_description.hpp"
#include "NLChain_4_problem_description.hpp"
#include "NLChain_6_problem_description.hpp"
#include "NLChain_8_problem_description.hpp"
#include "NLChain_10_problem_description.hpp"
#include "NLChain_12_problem_description.hpp"
#include "NLChain_14_problem_description.hpp"

using namespace grampc;
using namespace std::placeholders;

int main()
{
  typeInt numberOfChainElements = 6;
  typeInt Nx = 6 * numberOfChainElements - 3;

  // State  and parameter distribution
  Vector x0 = Vector::Zero(Nx);
  for (typeInt i = 0; i < numberOfChainElements; ++i)
  {
    x0[i * 6] = (i + 1) * 0.55;
  }
  Matrix P0 = 1e-6 * Matrix::Identity(Nx, Nx);
  DistributionPtr state = Gaussian(x0, P0);

  // Point transformation that can be used for the uncertainty propagation
  PointTransformationPtr transform = UT(state->dimension(), state->dimension(), 1, 2, 0);

  // Prediction horizon, sampling time
  typeRNum Thor = 3;
  typeRNum dt_MPC = 0.1;
  typeRNum dt_simulation = 0.01;
  typeRNum Tsim = 12;

  // Desired values for states and inputs
  std::vector<typeRNum> xdes(Nx, 0.0);
  std::vector<typeRNum> udes = {0, 0, 0};

  // Input constraints
  std::vector<typeRNum> umax = {10.0, 10.0, 10.0};
  std::vector<typeRNum> umin = {-10.0, -10.0, -10.0};

  // Declaration of vectors for the problem description
  std::vector<typeRNum> pCost = {10, 1, 0.01};

  // Initial input and tolerance of the constraints
  std::vector<typeRNum> u0 = {0, 0, 0};

  // create nominal problem description
  ProblemDescriptionPtr chainProblem;
  if (numberOfChainElements == 2)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_2_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 3)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_3_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 4)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_4_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 6)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_6_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 8)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_8_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 10)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_10_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 12)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_12_ProblemDescription(pCost));
  }
  else if (numberOfChainElements == 14)
  {
    chainProblem = ProblemDescriptionPtr(new Chain_14_ProblemDescription(pCost));
  }
  else
  {
    std::cerr << "There is no problem description for this number of chain elements!" << std::endl;
  }

  // configure stochastic problem description
  SigmaPointProblemDescriptionPtr problem = SigmaPointProblem(chainProblem, transform);

  // create solver
  GrampcPtr solver = Solver(problem);
  const typeGRAMPCparam *par = solver->getParameters();

  // simulator
  SystemFct trueSystemFunction = std::bind(&ProblemDescription::ffct, chainProblem, _1, _2, _3, _4, _5, par);
  Simulator sim(state->mean(), (typeInt) u0.size(), trueSystemFunction, "heun", 0, dt_MPC, dt_simulation, true);

  // set initial states and parameters depending on the propagation method
  problem->compute_x0_and_p0(state);
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
  solver->setopt_int("MaxGradIter", 3);
  solver->setopt_int("MaxMultIter", 1);
  solver->setopt_int("Nhor", 80);
  solver->setopt_string("Integrator", "erk2");
  solver->setopt_string("TerminalCost", "off");

  // computation time
  typeRNum compTime_ms = 0.0;

  // MPC loop
  for (typeRNum t = 0; t < Tsim; t += dt_MPC)
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
  std::cout << "Average computation time per sampling step: " << compTime_ms / std::round(Tsim / dt_MPC) << " ms" << std::endl;
  std::cout << "Program finished." << std::endl;
  return 0;
}