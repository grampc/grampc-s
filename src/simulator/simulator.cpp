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


#include "simulator/simulator.hpp"

namespace grampc
{
    Simulator::Simulator(VectorConstRef initialState, VectorConstRef param, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator, MatrixConstRef diffMatrixWienerProcess,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile)
    : systemFct_(systemFct),
      integrator_(integrator),
      state_(initialState),
      nextPoint_(initialState),
      param_(param),
      t_(t0),
      dt_MPC_(dt_MPC),
      dt_simulation_(dt_simulation),
      numStates_((typeInt) initialState.rows()),
      numParams_((typeInt) param.rows()),
      numControls_(numberOfControls),
      timeDerivative_(numStates_),
      timeDerivative2_(numStates_),
      diffMatrixWienerProcess_(diffMatrixWienerProcess),
      writeDataToFile_(writeDataToFile),
      noise_(numStates_)
    {
        // Check user input
        if(dt_MPC < dt_simulation)
        {
            std::cerr << "Sampling time of the simulation is smaller than sampling time of the controller." << std::endl;
        }
        else if(integrator_.compare("euler") && integrator_.compare("Euler") && integrator_.compare("heun") && integrator_.compare("Heun")
                && integrator_.compare("Euler-Maruyama") && integrator_.compare("euler-maruyama"))
        {
            std::cerr << "Integrator is unknown." << std::endl;
        }
        else if(!integrator_.compare("Euler"))
        {
            integrator_ = "euler";
        }
        else if(!integrator_.compare("Heun"))
        {
            integrator_ = "heun";
        }
        else if(!integrator_.compare("Euler-Maruyama"))
        {
            integrator_ = "euler-maruyama";
        }
        else if(std::fmod(dt_MPC, dt_simulation) > 1e-2)
        {
            std::cerr << "The sampling time of the model predictive controller must be a multiple of the simulation sampling time." << std::endl;
        }

        // write vector dimensions to file
        if(writeDataToFile_)
        {
            this->writeDimensionData();
        }

        // For stochastic simulations
        rng_ = RandomNumberGenerator(0);
    }

    Simulator::Simulator(VectorConstRef initialState, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator, MatrixConstRef diffMatrixWienerProcess,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile)
    : Simulator(initialState, Vector::Zero(0), numberOfControls, systemFct, integrator, diffMatrixWienerProcess, t0, dt_MPC, dt_simulation, writeDataToFile)
    {
    }

    Simulator::Simulator(VectorConstRef initialState, VectorConstRef param, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile)
    : Simulator(initialState, param, numberOfControls, systemFct, integrator, Matrix::Zero(0,0), t0, dt_MPC, dt_simulation, writeDataToFile)
    {
    }

    Simulator::Simulator(VectorConstRef initialState, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile)
    : Simulator(initialState, Vector::Zero(0), numberOfControls, systemFct, integrator, Matrix::Zero(0,0), t0, dt_MPC, dt_simulation, writeDataToFile)
    {
    }

    const Vector &Simulator::simulate(VectorConstRef control)
    {
        // end of simulation
        typeRNum simulationEndTime = t_ + dt_MPC_;

        if (integrator_ == "euler")
        {
            while (t_ < simulationEndTime)
            {
                // evaluate system function
                systemFct_(timeDerivative_, t_, state_, control, param_);

                // x_{k+1} = x_k + dt*x'_k
                state_ += dt_simulation_ * timeDerivative_;

                // update time
                t_ += dt_simulation_;
            }
        }
        else if (integrator_ == "heun")
        {
            while (t_ < simulationEndTime)
            {
                // evaluate system function
                systemFct_(timeDerivative_, t_, state_, control, param_);

                // compute Euler: x_{k+1} = x_k + dt*x'_k
                nextPoint_ = state_ + dt_simulation_ * timeDerivative_;

                // evaluate system function at the second point
                systemFct_(timeDerivative2_, t_ + dt_simulation_, nextPoint_, control, param_);

                // compute Heun: x_{k+1} = x_k + 0.5 * dt * (x1' + x2')
                state_ += 0.5 * dt_simulation_ * (timeDerivative_ + timeDerivative2_);

                // update time
                t_ += dt_simulation_;
            }
        }
        else if (integrator_ == "euler-maruyama")
        {
             while (t_ < simulationEndTime)
            {
                // evaluate system function
                systemFct_(timeDerivative_, t_, state_, control, param_);

                // sample from Gaussian distribution
                for(typeInt i = 0; i < numStates_; ++i)
                {
                    noise_(i) = dist_(rng_);
                }

                // x_{k+1} = x_k + dt*x'_k + Sigma*W_k
                state_ += dt_simulation_ * timeDerivative_ + diffMatrixWienerProcess_ * std::sqrt(dt_simulation_) * noise_;

                // update time
                t_ += dt_simulation_;
            }
        }

        // write data to file
        if(writeDataToFile_)
        {
            this->writeClosedLoopData(control);
        }

        // return state
        return state_;
    }

    typeRNum Simulator::time()
    {
        return t_;
    }

    void Simulator::writeDimensionData()
    {
        // open output files and remove data of privious simulations
        std::ofstream tout("tvec.txt");
        std::ofstream xout("xvec.txt");
        std::ofstream uout("uvec.txt");
        std::ofstream dimOut("dim.txt");

        // write dimensions
        dimOut << "VariableName,Data" << std::endl;
        dimOut << "Nx," << numStates_ << std::endl;
        dimOut << "Nu," << numControls_ << std::endl;

        // time
        tout << t_  << std::endl;

        // states
        for(typeInt i = 0; i < numStates_; ++i)
        {
            xout << state_.row(i) << "\t";
        }
        xout << std::endl;

        // close file
        dimOut.close();
        uout.close();
        xout.close();
        tout.close();
    }

    void Simulator::writeClosedLoopData(VectorConstRef control)
    {
        // open output files
        std::ofstream tout("tvec.txt", std::ios::app);
        std::ofstream xout("xvec.txt", std::ios::app);
        std::ofstream uout("uvec.txt", std::ios::app);

        // time
        tout << t_  << std::endl;

        // states
        for(typeInt i = 0; i < numStates_; ++i)
        {
            xout << state_.row(i) << "\t";
        }
        xout << std::endl;

        // controls
        for(typeInt i = 0; i < numControls_; ++i)
        {
            uout << control[i] << "\t";
        }
        uout << std::endl;

        // close files
        uout.close();
        xout.close();
        tout.close();
    }
}