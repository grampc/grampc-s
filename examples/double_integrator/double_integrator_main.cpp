#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include "grampc_s.hpp"
#include "double_integrator_problem_description.hpp"

using namespace grampc;

int main()
{
    // State  and parameter distribution
    Vector x0(2);
    x0 << 2, 0;
    DistributionPtr state = Gaussian(x0, 1e-6 * Matrix::Identity(2, 2));

    Vector WienerProcessVariance(2);
    WienerProcessVariance << 0.001, 0.001;
    Matrix covProcessNoise = WienerProcessVariance.asDiagonal();

    // Point transformation that can be used for the uncertainty propagation
    PointTransformationPtr transform = UT(state->dimension(), state->dimension(), 0.1, 2, 0);
    // PointTransformationPtr transform = StirlingFirstOrder(state->dimension(), state->dimension(), 1);
    // PointTransformationPtr transform = StirlingSecondOrder(state->dimension(), state->dimension(), 1);

    // Prediction horizon, sampling time, and initial time
    typeRNum Thor = 4;
    typeRNum dt = 0.001;

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

    // create nominal problem description
    StochasticProblemDescriptionPtr integratorProblem = StochasticProblemDescriptionPtr(new DoubleIntegratorProblemDescription(pCost, pCon));

    // configure stochastic problem description
    grampc::ResamplingProblemDescription problem(integratorProblem, constraintApprox, transform, covProcessNoise);

    // create solver
    grampc::Grampc solver(&problem);
    const typeGRAMPCparam *par = solver.getParameters();

    // set initial states and parameters depending on the propagation method
    problem.compute_x0_and_p0(state);
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
    solver.setopt_int("MaxMultIter", 3);
    solver.setopt_int("Nhor", 20);
    solver.setopt_string("Integrator", "heun");
    solver.setopt_real("PenaltyMin", 1e3);
    solver.setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
    solver.setopt_string("OptimTime", "on");

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
    dimOut << "Np," << 0 << std::endl;
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
