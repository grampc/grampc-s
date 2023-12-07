#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "grampc_s.hpp"
#include "reactor_problem_description.hpp"

using namespace grampc;

int main()
{
     // Parameter for the propagation methods
    typeInt MC_numPoints = 100;
    RandomNumberGenerator rng = RandomNumberGenerator(1);
    typeRNum UT_alpha = 0.1;
    typeRNum UT_beta = 2.0;
    typeRNum UT_kappa = 0;

    // State  and parameter distribution
    Vector x0(2);
    x0 << 0.708, 0.09;
    DistributionPtr state = Gaussian(x0, 1e-3 * Matrix::Identity(2, 2));
    MultiDistributionPtr param = MultiDist({Gaussian(0.0, 1e-6)});

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = UT(state->dimension() + param->dimension(), UT_alpha, UT_beta, UT_kappa, {1, 1, 0});
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension() + param->dimension(), 1, {1, 1, 0});
    // PointTransformationPtr transform = StirlingSecondOrder(state->dimension() + param->dimension(), 1, {1, 1, 0});
    // PointTransformationPtr transform = MonteCarlo(state->dimension() + param->dimension(), MC_numPoints, rng);
    // PointTransformationPtr transform = Quadrature(MultiDist({state, param})->polynomialFamily(), {1, 1, 0});
    // PointTransformationPtr transform = PCE(MultiDist({state, param})->polynomialFamily(), 2, {1, 1, 0});

    // Prediction horizon, sampling time, and initial time
    typeRNum Thor = 0.1;
    typeRNum dt = 0.01;

    // Desired values for states and inputs
    std::vector<typeRNum> xdes = {0.215, 0.09};
    std::vector<typeRNum> udes = {19.6};

    // Input constraints
    std::vector<typeRNum> umax = {400};
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
    pCost.push_back(1.0); // Coefficient for (x1-x1_Des)^2
    pCost.push_back(1.0); // Coefficient for (x2-x2_Des)^2

    // Input cost R
    pCost.push_back(0.1); // Coefficient for (u-u_Des)^2

    // Constraint parameters
    pCon.push_back(0.14);

    // Chance constraint probability
    Vector vec_chance_constraint(1);
    vec_chance_constraint << 0.9;

    // Chance constraint approximation
    // ChebyshevConstraintApproximation constraintApprox(&vec_chance_constraint);
    ChanceConstraintApproximationPtr constraintApprox = GaussianApprox(vec_chance_constraint);
    // SymmetricConstraintApproximation constraintApprox(&vec_chance_constraint);

    // create nominal problem description
    StochasticProblemDescriptionPtr reactor_problem = StochasticProblemDescriptionPtr(new ReactorProblemDescription(pSys, pCost, pCon));

    // configure stochastic problem description
    grampc::SigmaPointProblemDescription problem(reactor_problem, constraintApprox, transform);
    // grampc::MonteCarloProblemDescription problem(reactor_problem, constraintApprox, transform);

    // create solver
    grampc::Grampc solver(&problem);
    const typeGRAMPCparam *par = solver.getParameters();

    // set initial states and parameters depending on the propagation method
    problem.compute_x0_and_p0(state, param);
    solver.setparam_real_vector("x0", problem.x0());
    solver.setparam_real_vector("p0", problem.p0());

    xdes.resize(par->Nx, 0.0);
    constraintsAbsTol.resize(par->Nc, 0.0);

    // set parameters
    solver.setparam_real("Thor", Thor);
    solver.setparam_real("dt", dt);
    solver.setparam_real_vector("xdes", &xdes[0]);
    solver.setparam_real_vector("u0", &u0[0]);
    solver.setparam_real_vector("udes", &udes[0]);
    solver.setparam_real_vector("umax", &umax[0]);
    solver.setparam_real_vector("umin", &umin[0]);

    // set options
    solver.setopt_int("MaxGradIter", 3);
    solver.setopt_int("MaxMultIter", 2);
    solver.setopt_int("Nhor", 20);
    solver.setopt_string("Integrator", "heun");
    solver.setopt_real("PenaltyMin", 3e7);
    solver.setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
    solver.setopt_string("InequalityConstraints", "on");
    solver.setopt_string("EqualityConstraints", "off");
    solver.setopt_string("ConvergenceCheck", "off");
    solver.setopt_string("ShiftControl", "off");

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
    typeInt Nx = solver.getParameters()->Nx;
    typeInt Nu = solver.getParameters()->Nu;
    typeInt Nc = solver.getParameters()->Nc;

    // run gradient iterations and save computation time
    std::chrono::steady_clock::time_point tic = std::chrono::steady_clock::now();
    solver.run();
    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();
    compTime_ms += std::chrono::duration_cast<std::chrono::microseconds>(toc - tic).count() / 1000.0;

    // write dimensions
    dimOut << "VariableName,Data" << std::endl;
    dimOut << "Nx," << state->dimension() << std::endl;
    dimOut << "Nu," << u0.size() << std::endl;
    dimOut << "Np," << param->dimension() << std::endl;
    dimOut << "Nh," << solver.getParameters()->Nh;

    typeRNum tempConstraint[Nc];

    for (typeInt i = 0; i < solver.getOptions()->Nhor; ++i)
    {
        tout << solver.getWorkspace()->t[i] << std::endl;
        for (typeInt j = 0; j < Nx; ++j)
        {
            xout << solver.getWorkspace()->x[i * Nx + j] << "\t";
            adjout << solver.getWorkspace()->adj[i * Nx + j] << "\t";
        }
        xout << std::endl;
        adjout << std::endl;
        for (typeInt j = 0; j < Nu; ++j)
        {
            uout << solver.getWorkspace()->u[i * Nu + j] << "\t";
        }
        uout << std::endl;

        problem.hfct(tempConstraint, solver.getWorkspace()->t[i], solver.getWorkspace()->x + i * Nx, solver.getWorkspace()->u + i * Nu, solver.getWorkspace()->p);
        for (typeInt j = 0; j < Nc; ++j)
        {
            conOut << tempConstraint[j] << "\t";
        }
        conOut << std::endl;
    }

    // Print computation time
    std::cout << "Computation time: " << compTime_ms << " ms" << std::endl;

    std::cout << "Program finished" << std::endl;
    return 0;
}
