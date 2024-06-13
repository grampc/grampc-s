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


#ifndef GAUSSIAN_PROCESS_HPP
#define GAUSSIAN_PROCESS_HPP

#include <iostream>
#include <fstream>
#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    // Zero-mean Gaussian processes with univariate output
    class GaussianProcess
    {
    public:
        // Constructor of a Gaussian process with a fixed number of data points 
        GaussianProcess(const GaussianProcessData& data, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency);

        // Constructor of a Gaussian process using data from a text file
        GaussianProcess(const std::string& fileName, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency);

        // Mean of the Gaussian process at the evaluation point
        typeRNum mean(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl);

        // Variance of the Gaussian process at the evaluation point
        typeRNum variance(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl);

        // Gradient dmean/dx multiplied by scalar vec, i.e. (dmean/dx)^T*vec
        const Vector& dmean_dx_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec);

        // Gradient dmean/du multiplied by scalar vec, i.e. (dmean/du)^T*vec
        const Vector& dmean_du_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec);

        // Gradient dvar/dx multiplied by scalar vec, i.e. (dvar/dx)^T*vec
        const Vector& dvar_dx_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec);

        // Gradient dvar/du multiplied by scalar vec, i.e. (dvar/du)^T*vec
        const Vector& dvar_du_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec);

    private:  
        // Read input and output data for a Gaussian process from a text file
        GaussianProcessData readGaussianProcessData(const std::string& fileName) const;

        // Compute difference between an evaluation point and all data points
        void pointDifference(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl);

        typeInt inputDim_;
        typeInt numDataPoints_;
        const std::vector<bool> stateDependency_;
        const std::vector<bool> controlDependency_;
        typeInt numStates_;
        typeInt numControls_;
        typeInt numStateArguments_;
        typeInt numControlArguments_;
        Matrix inputData_;
        StationaryKernelConstPtr kernel_;
        Matrix K_inv_;
        Vector K_inv_y_;
        Vector k_evaluationPoint_;
        Matrix pointDiffMatrix_;
        Vector temp_vec_numDataPoints_;
        Matrix dk_dx_;
        Matrix dk_du_;
        Vector dmean_dx_vec_;
        Vector dmean_du_vec_;
        Vector dvar_dx_vec_;
        Vector dvar_du_vec_;
        IntVector indicesStates_;
        IntVector indicesControls_;
        Vector evaluationPoint_;
        typeRNum kernelZero_;
    };

    // Alias
    typedef std::shared_ptr<GaussianProcess> GaussianProcessPtr;
    typedef std::shared_ptr<const GaussianProcess> GaussianProcessConstPtr;

    // Constructor of a Gaussian process with a fixed number of data points 
    GaussianProcessPtr GP(const GaussianProcessData& data, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency);
    // Constructor of a Gaussian process using data from a text file
    GaussianProcessPtr GP(const std::string& fileName, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency);
}

#endif // GAUSSIAN_PROCESS_HPP
