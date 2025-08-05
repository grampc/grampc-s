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


#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "util/grampc_s_constants.hpp"
#include "point_transformation/monte_carlo.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    class Simulator 
    {
    public:
        // Constructor for a simulator with states and parameters 
        Simulator(VectorConstRef initialState, VectorConstRef param, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator, MatrixConstRef diffMatrixWienerProcess,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile);

        // Constructor for a simulator with states only
        Simulator(VectorConstRef initialState, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator, MatrixConstRef diffMatrixWienerProcess,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile);

        // Constructor for a simulator with states and parameters and without a Wiener process
        Simulator(VectorConstRef initialState, VectorConstRef param, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile);

        // Constructor for a simulator with states only and without a Wiener process
        Simulator(VectorConstRef initialState, typeInt numberOfControls, const SystemFct &systemFct, const std::string& integrator,
                  typeRNum t0, typeRNum dt_MPC, typeRNum dt_simulation, typeBoolean writeDataToFile);

        // Simulate system and return state
        const Vector& simulate(VectorConstRef control);

        // Return current time
        typeRNum time();

    private:
        // write vector dimensions to file
        void writeDimensionData();

        // write time, state and control data to files
        void writeClosedLoopData(VectorConstRef control);

        SystemFct systemFct_;
        typeRNum t_;
        typeRNum dt_MPC_;
        typeRNum dt_simulation_;
        typeInt numStates_;
        typeInt numParams_;
        typeInt numControls_;
        std::string integrator_;
        Vector state_;
        Vector nextPoint_;
        Vector param_;
        Vector timeDerivative_;
        Vector timeDerivative2_;
        Matrix diffMatrixWienerProcess_;
        bool writeDataToFile_;
        std::normal_distribution<typeRNum> dist_{0.0, 1.0};
        RandomNumberGenerator rng_;
        Vector noise_;
    };
}


#endif // SIMULATOR_HPP