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


#include "gaussian_process/gaussian_process.hpp"

namespace grampc
{
    GaussianProcess::GaussianProcess(const GaussianProcessData& data, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency)
    : inputDim_(data.inputDimension),
      numDataPoints_(data.numberOfDataPoints),
      stateDependency_(stateDependency),
      controlDependency_(controlDependency),
      numStates_((typeInt) stateDependency.size()),
      numControls_((typeInt) controlDependency.size()),
      numStateArguments_((typeInt) std::count(stateDependency.begin(), stateDependency.end(), true)),
      numControlArguments_((typeInt) std::count(controlDependency.begin(), controlDependency.end(), true)),
      inputData_(data.inputData),
      kernel_(kernel),
      K_inv_(numDataPoints_, numDataPoints_),
      K_inv_y_(numDataPoints_),
      k_evaluationPoint_(numDataPoints_),
      pointDiffMatrix_(inputDim_, numDataPoints_),
      temp_vec_numDataPoints_(numDataPoints_),
      dk_dx_(numStateArguments_, numDataPoints_),
      dk_du_(numControlArguments_, numDataPoints_),
      dmean_dx_vec_(Vector::Zero(numStates_)),
      dmean_du_vec_(Vector::Zero(numControls_)),
      dvar_dx_vec_(Vector::Zero(numStates_)),
      dvar_du_vec_(Vector::Zero(numControls_)),
      indicesStates_(numStateArguments_),
      indicesControls_(numControlArguments_),
      evaluationPoint_(inputDim_),
      kernelZero_(kernel->evaluate(Vector::Zero(inputDim_)))
    {
        // Check user input
        if(inputData_.cols() != data.outputData.rows())
        {
            std::cerr << "The number of input data points is not equal to the number of output data points!" << std::endl;
        }else if(inputDim_ != kernel->inputDimension())
        {
            std::cerr << "The input dimension of the kernel does not match the dimension of the input data!" << std::endl;
        }else if(numControlArguments_ + numStateArguments_ != inputDim_)
        {
            std::cerr << "The number of rows of the input data is not equal to the number of state and control arguments of the GP!" << std::endl;
        }
        
        // covariance matrix and temporary vector tau
        Matrix K(numDataPoints_, numDataPoints_);
        Vector tau(inputDim_);

        // Compute the elements of the covariance matrix
        for(typeInt i = 0; i < numDataPoints_;  ++i)
        {
            for(typeInt j = 0; j < numDataPoints_; ++j)
            {
                tau = inputData_.col(i) - inputData_.col(j);
                K(i, j) = kernel_->evaluate(tau);
            }
        }

        // Add output noise to the diagonal elements
        K += data.outputNoiseVariance * Matrix::Identity(numDataPoints_, numDataPoints_);
        
        // QR decomposition of the covariance matrix, check that the inverse exist
        Eigen::ColPivHouseholderQR<Matrix> dec(K);
        if(!dec.isInvertible())
        {
            std::cerr << "Matrix K is not invertible!" << std::endl;
        }

        // Compute inverse matrix and product of inverse and outputdata 
        K_inv_  = dec.inverse();
        K_inv_y_ = K_inv_ * data.outputData;

        // Indices of states and controls in the input data
        for(typeInt i = 0; i < numStateArguments_; ++i)
        {
            indicesStates_(i) = i;
        }
        for(typeInt i = numStateArguments_; i < inputDim_; ++i)
        {
            indicesControls_(i - numStateArguments_) = i;
        }
    }

    GaussianProcess::GaussianProcess(const std::string& fileName, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency)
    : GaussianProcess(readGaussianProcessData(fileName), kernel, stateDependency, controlDependency)
    {
    }

    typeRNum GaussianProcess::mean(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            k_evaluationPoint_(i) = kernel_->evaluate(pointDiffMatrix_.col(i));
        }
        return k_evaluationPoint_.transpose() * K_inv_y_;
    }

    typeRNum GaussianProcess::variance(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            k_evaluationPoint_(i) = kernel_->evaluate(pointDiffMatrix_.col(i));
        }
        temp_vec_numDataPoints_.noalias() = K_inv_ * k_evaluationPoint_;
        return kernelZero_ - k_evaluationPoint_.transpose() * temp_vec_numDataPoints_;
    }

    const Vector& GaussianProcess::dmean_dx_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute transposed Jacobian of vector k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            kernel_->gradient(dk_dx_.col(i), pointDiffMatrix_.col(i), indicesStates_);
        }

        // gradient of the mean at the evaluation point if GP depends on this state, else return 0
        typeInt index = 0;
        for(typeInt i = 0; i < numStates_; ++i)
        {
            if(stateDependency_[i])
            {
                dmean_dx_vec_(i) = (dk_dx_.row(index++) * K_inv_y_).value() * vec;
            }
        }
        return dmean_dx_vec_;
    }

    const Vector& GaussianProcess::dmean_du_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute transposed Jacobian of vector k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            kernel_->gradient(dk_du_.col(i), pointDiffMatrix_.col(i), indicesControls_);
        }

        // gradient of the mean at the evaluation point if GP depends on this input, else return 0
        typeInt index = 0;
        for(typeInt i = 0; i < numControls_; ++i)
        {
            if(controlDependency_[i])
            {
                dmean_du_vec_(i) = (dk_du_.row(index++) * K_inv_y_).value() * vec;
            }
        }
        return dmean_du_vec_;
    }

    const Vector& GaussianProcess::dvar_dx_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute k(x*, X) and transposed Jacobian of vector k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            k_evaluationPoint_(i) = kernel_->evaluate(pointDiffMatrix_.col(i));
            kernel_->gradient(dk_dx_.col(i), pointDiffMatrix_.col(i), indicesStates_);
        }

        // return -2 * dk_dx^T(x*, X) * K(X, X)^(-1) * k(x*, X) * vec if GP depends on this state, else return 0
        typeInt index = 0;
        for(typeInt i = 0; i < numStates_; ++i)
        {
            if(stateDependency_[i])
            {
                temp_vec_numDataPoints_.noalias() = K_inv_ * k_evaluationPoint_;
                dvar_dx_vec_(i) = -2.0 * (dk_dx_.row(index++) * temp_vec_numDataPoints_).value() * vec;
            }
        }
        return dvar_dx_vec_;
    }

    const Vector& GaussianProcess::dvar_du_vec(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl, typeRNum vec)
    {
        // difference between the evaluation point and all data points
        pointDifference(evaluationPointState, evaluationPointControl);

        // Compute k(x*, X) and transposed Jacobian of vector k(x*, X)
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            k_evaluationPoint_(i) = kernel_->evaluate(pointDiffMatrix_.col(i));
            kernel_->gradient(dk_du_.col(i), pointDiffMatrix_.col(i), indicesControls_);
        }

        // return -2 * dk_dx^T(x*, X) * A * k(x*, X) if GP depends on this input, else return 0
        typeInt index = 0;
        for(typeInt i = 0; i < numControls_; ++i)
        {
            if(controlDependency_[i])
            {
                temp_vec_numDataPoints_.noalias() = K_inv_ * k_evaluationPoint_;
                dvar_du_vec_(i) = -2.0 * (dk_du_.row(index++) * temp_vec_numDataPoints_).value() * vec;
            }
        }
        return dvar_du_vec_;
    }

    GaussianProcessData GaussianProcess::readGaussianProcessData(const std::string& fileName) const
    {
        GaussianProcessData data;
        
        // define input stream and read matrix dimensions and output noise variance
        std::ifstream inputStream(fileName);
        if(inputStream.fail()) 
        {
            std::cerr << "Failed to open GP-data file." <<std::endl;
        }
        inputStream >> data.inputDimension >> data.numberOfDataPoints >> data.outputNoiseVariance;

        // resize matrices with dimensions from the text file
        data.inputData.resize(data.inputDimension, data.numberOfDataPoints);
        data.outputData.resize(data.numberOfDataPoints);

        // read input data for the Gaussian process
        for (typeInt row = 0; row < data.inputDimension; ++row)
        {
            for (typeInt col = 0; col < data.numberOfDataPoints; ++col)
            {
                inputStream >> data.inputData(row, col);
            }
        }

        // read output data for the Gaussian process
        for (typeInt col = 0; col < data.numberOfDataPoints; ++col)
        {
            inputStream >> data.outputData(col);
        }

        // close file and return data
        inputStream.close();
        return data;
    }

    void GaussianProcess::pointDifference(VectorConstRef evaluationPointState, VectorConstRef evaluationPointControl)
    {
        // Compute evaluation point in correct data format
        typeInt index = 0;
        for(typeInt i = 0; i < numStates_; ++i)
        {
            if(stateDependency_[i])
            {
                evaluationPoint_(index++) = evaluationPointState(i);
            }
        }
        for(typeInt i = 0; i < numControls_; ++i)
        {
            if(controlDependency_[i])
            {
                evaluationPoint_(index++) = evaluationPointControl(i);
            }
        }

        // difference between the evaluation point and all data points
        for(typeInt i = 0; i < numDataPoints_; ++i)
        {
            pointDiffMatrix_.col(i) = evaluationPoint_ - inputData_.col(i);
        }
    }

    GaussianProcessPtr GP(const GaussianProcessData& data, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency)
    {
        return GaussianProcessPtr(new GaussianProcess(data, kernel, stateDependency, controlDependency));
    }

    GaussianProcessPtr GP(const std::string& fileName, StationaryKernelConstPtr kernel, const std::vector<bool>& stateDependency, const std::vector<bool>& controlDependency)
    {
        return GaussianProcessPtr(new GaussianProcess(fileName, kernel, stateDependency, controlDependency));
    }
}