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


#include "problem_description/resampling_GP_problem_description.hpp"

namespace grampc
{
    ResamplingGPProblemDescription::ResamplingGPProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                                    PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP)
        : numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation),
          gaussianProcessVec_(gaussianProcessVec),
          dynamicsIndicesWithGP_(dynamicsIndicesWithGP),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient()),
          tempScalar_(1)
        {
        typeInt Ng, NgT;

        // Call ocp_dim to get the number of states, controls, parameters, and constraints
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, &Ng, &numConstraints_, &NgT, &numTerminalConstraints_); 

        // Dimension of each point
        pointDim_ = numStates_ + numParams_;

        // Distribution of states and parameters
        meanStateAndParam_ = Vector::Zero(pointDim_);
        covStateAndParam_ = Matrix::Zero(pointDim_, pointDim_);

        // Allocate memory for the Cholesky decomposition
        llt_.compute(covStateAndParam_);

        // allocate memory
        x0_ = Vector::Zero(numStates_ + numStates_ * pointDim_);
        p0_ = Vector::Zero(numParams_);
        constraintStdDev_ = Vector::Zero(numConstraints_);
        terminalConstraintStdDev_ = Vector::Zero(numTerminalConstraints_);
        pointsTransformed_ = Matrix::Zero(numStates_, numSigmaPoints_);
        covStateAndParam_ = Matrix::Zero(pointDim_, pointDim_);
        dcov_dPointsY_ = Vector::Zero(numStates_ * numSigmaPoints_);
        constraintMatrix_ = Matrix::Zero(numConstraints_, numSigmaPoints_);
        terminalConstraintMatrix_ = Matrix::Zero(numTerminalConstraints_, numSigmaPoints_);
        temp_vec_pointDim_numPoints_ = Vector::Zero(numSigmaPoints_  *pointDim_);
        temp_vec_numInputs_ = Vector::Zero(numControlInputs_);
        dmean1D_dPoints_ = pointTransformation_->dmean1D_dpoints().transpose();
        constraintVec_ = Matrix::Zero(numConstraints_, numSigmaPoints_);
        terminalConstraintVec_ = Matrix::Zero(numTerminalConstraints_, numSigmaPoints_);

        // Array dimensions
        Nx_ = numStates_ + numStates_ * numStates_ + numStates_ * numParams_;
        Np_ = numParams_;
        Nu_ = numControlInputs_;
        Nh_ = numConstraints_;
        NhT_ = numTerminalConstraints_;
    }

    ResamplingGPProblemDescription::ResamplingGPProblemDescription(ProblemDescriptionPtr problemDescription, 
                                                               PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP)
    : ResamplingGPProblemDescription(problemDescription, Chebyshev(Vector(0)), pointTransformation, gaussianProcessVec, dynamicsIndicesWithGP)
    {
        if(numConstraints_ > 0)
        {
            std::cerr << "Constructor without chance constraint approximation is called but constraints exist in the problem description!" << std::endl;
        }
    }

    void ResamplingGPProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Np = Np_;
        *Nu = Nu_;
        *Nh = Nh_;
        *NhT = NhT_;
        *Ng = 0;
        *NgT = 0;
    }

    void ResamplingGPProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Mapping of the covariance output
        Eigen::Map<Matrix> d_cov(out.data() + numStates_, pointDim_, numStates_);

        // Mapping of the covariance matrix
        Eigen::Map<const Matrix> cov(x.data() + numStates_, pointDim_, numStates_);

        typeInt index;
        
        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = x.segment(0, numStates_);
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Transform points
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Known dynamics
            problemDescription_->ffct(pointsTransformed_.col(i), t, pointsInitial_.col(i), u, pointsInitial_.col(i).segment(numStates_, numParams_), param);

            // GP-dynamics
            index = 0;
            for(typeInt j : dynamicsIndicesWithGP_)
            {
                pointsTransformed_(j,i) += gaussianProcessVec_[index++]->mean(pointsInitial_.col(i), u);
            }
        }

        // Cross-covariance between initial points and transformed points
        const Matrix& crossCov = pointTransformation_->covariance(pointsInitial_, pointsTransformed_);

        // Time derivative of the mean
        out.segment(0, numStates_) = pointTransformation_->mean(pointsTransformed_);

        // Time derivative of the covariance
        d_cov = crossCov;
        d_cov.topRows(numStates_) += crossCov.topRows(numStates_).transpose();

        // variance from the Gaussian process
        index = 0;
        for (typeInt i : dynamicsIndicesWithGP_)
        {
            d_cov(i, i) += gaussianProcessVec_[index++]->variance(pointsInitial_.col(0), u);
        }
    }

    void ResamplingGPProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param)
    {   
        typeInt index;

        // Mapping of the covariance matrix
        Eigen::Map<const Matrix> cov(x.data() + numStates_, pointDim_, numStates_);

        // Mapping of the covariance output
        Eigen::Map<Matrix> dd_cov(out.data() + numStates_, pointDim_, numStates_);

         // Mapping of adj
        Eigen::Map<const Vector> adjCov(adj.data() + numStates_, pointDim_ * numStates_);
        
        // covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = x.segment(0, numStates_);
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Transform points
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Known dynamics
            problemDescription_->ffct(pointsTransformed_.col(i), t, pointsInitial_.col(i), u, pointsInitial_.col(i).segment(numStates_, numParams_), param);

            // GP-dynamics
            index = 0;
            for(typeInt j : dynamicsIndicesWithGP_)
            {
                pointsTransformed_(j,i) += gaussianProcessVec_[index++]->mean(pointsInitial_.col(i), u);
            }
        }

        // d(d_mean_state)/d(mean_state)
        const Vector& dmean_dpoints_vec = pointTransformation_->dmean_dpoints_vec(adj.segment(0, numStates_));
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_, numStates_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_),
                dmean_dpoints_vec.segment(i*numStates_, numStates_),
                param
            );
            problemDescription_->dfdp_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_ + numStates_, numParams_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                dmean_dpoints_vec.segment(i*numStates_, numStates_),
                param
            );
        
            // GP-dynamics
            index = 0;
            for(typeInt j : dynamicsIndicesWithGP_)
            {
                temp_vec_pointDim_numPoints_.segment(i * pointDim_, numStates_) += gaussianProcessVec_[index++]->dmean_dx_vec(pointsInitial_.col(i), u, dmean_dpoints_vec(i * numStates_ + j));
            }
        }
        out.segment(0, numStates_) = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);

        // d(d_mean_state)/d(cov)
        dd_cov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);

        // d(d_cov)/d(mean)
        dcov_dPointsY_ = pointTransformation_->dcov_dpointsY_vec(pointsInitial_, adjCov);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_, numStates_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_),
                dcov_dPointsY_.segment(i*numStates_, numStates_), 
                param
            );
            problemDescription_->dfdp_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_ + numStates_, numParams_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_),
                dcov_dPointsY_.segment(i*numStates_, numStates_), 
                param
            );
        
            // GP-dynamics
            index = 0;
            for(typeInt j : dynamicsIndicesWithGP_)
            {
                temp_vec_pointDim_numPoints_.segment(i * pointDim_, numStates_) += gaussianProcessVec_[index++]->dmean_dx_vec(pointsInitial_.col(i), u, dcov_dPointsY_(i * numStates_ + j));
            }
        }
        temp_vec_pointDim_numPoints_ += pointTransformation_->dcov_dpointsX_vec(pointsTransformed_, adjCov);
        
        // Influence of GP-variance
        index = 0;
        for(typeInt j : dynamicsIndicesWithGP_)
        {
            temp_vec_pointDim_numPoints_.segment(0, numStates_) += gaussianProcessVec_[index++]->dvar_dx_vec(pointsInitial_.col(0), u, adjCov(j + j * pointDim_));
        } 
        out.segment(0, numStates_) += pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);

        //d(d_cov)/d(cov)
        dd_cov += pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);
    }

    void ResamplingGPProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param)
    {
        typeInt index;
        out.setZero();

        // Mapping of the covariance matrix
        Eigen::Map<const Matrix> cov(x.data() + numStates_, pointDim_, numStates_);
        
        // covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = x.segment(0, numStates_);
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // d(d_mean)/du
        const Vector& dmean_dpoints_vec = pointTransformation_->dmean_dpoints_vec(adj);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Known dynamics
            problemDescription_->dfdu_vec(
                temp_vec_numInputs_, 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                dmean_dpoints_vec.segment(i*numStates_, numStates_), 
                param
            );
            out += temp_vec_numInputs_;
            
            // GP-dynamics
            index = 0;
            for(typeInt j : dynamicsIndicesWithGP_)
            {
                out += gaussianProcessVec_[index++]->dmean_du_vec(pointsInitial_.col(i), u, dmean_dpoints_vec(j + i * numStates_));
            }
        }

        // d(d_cov)/du
        dcov_dPointsY_ = pointTransformation_->dcov_dpointsY_vec(pointsInitial_, adj.segment(numStates_, pointDim_*numStates_));
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(
                temp_vec_numInputs_, 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                dcov_dPointsY_.segment(i*numStates_, numStates_),
                param
            );
            out += temp_vec_numInputs_;
        }
        // GP-dynamics
        index = 0;
        for (typeInt j : dynamicsIndicesWithGP_)
        {
            out += gaussianProcessVec_[index++]->dvar_du_vec(pointsInitial_.col(0), u, dcov_dPointsY_(j));
        }
    }

    void ResamplingGPProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        problemDescription_->lfct(out, t, x, u, p, param);
    }

    void ResamplingGPProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        problemDescription_->dldx(out, t, x, u, p, param);
    }

    void ResamplingGPProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        problemDescription_->dldu(out, t, x, u, p, param);
    }

    void ResamplingGPProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        problemDescription_->Vfct(out, t, x, p, param);
    }

    void ResamplingGPProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        problemDescription_->dVdx(out, t, x, p, param);
    }

    void ResamplingGPProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
       problemDescription_->dVdT(out, t, x, p, param);
    }

    void ResamplingGPProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> cov(x.data() + numStates_, pointDim_, numStates_);

        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = x.segment(0, numStates_);
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(constraintMatrix_.col(i), t, pointsInitial_.col(i), u,  pointsInitial_.col(i).segment(numStates_, numParams_), param);
        }
            
        // Compute tightened constraints
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            // Compute the standard deviation of the constraints
            constraintStdDev_(i) = std::sqrt(pointTransformation_->variance(constraintMatrix_.row(i)));

            // Tightened constraint
            out[i] = pointTransformation_->mean1D(constraintMatrix_.row(i)) + constraintTighteningCoeff_(i) * constraintStdDev_(i);
        }
    }

    void ResamplingGPProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        // Mapping of the outputs
        Eigen::Map<Matrix> outCov(out.data() + numStates_, pointDim_, numStates_);

        // Get points and Cholesky decomposition of the covariance matrix that are both computed in hfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();
        const Matrix& covCholStateAndParam_ = llt_.matrixLLT();

        // Compute the derivative of the tightened constraint with respect to the constraints of the sigma points
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            constraintVec_.row(i) = (dmean1D_dPoints_ + pointTransformation_->dvar_dpoints(constraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(i) / 2.0 / std::max(stdDevMin_, constraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhdx_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_, numStates_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                constraintVec_.col(i),
                param
            );
            problemDescription_->dhdp_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_ + numStates_, numParams_), 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                constraintVec_.col(i),
                param
            );
        }
        out.segment(0, numStates_) = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);
        outCov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);
    }

    void ResamplingGPProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        out.setZero();

        // Get points that are computed in hfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();

        // constraintVec is computed and saved in dhdx_vec()
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint with respect to the input
            problemDescription_->dhdu_vec(
                temp_vec_numInputs_, 
                t, 
                pointsInitial_.col(i), 
                u, 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                constraintVec_.col(i),
                param
            );
            out += temp_vec_numInputs_;
        }
    }

    void ResamplingGPProblemDescription::hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> cov(x.data() + numStates_, pointDim_, numStates_);

        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = x.segment(0, numStates_);
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(terminalConstraintMatrix_.col(i), t, pointsInitial_.col(i), pointsInitial_.col(i).segment(numStates_, numParams_), param);
        }
            
        // Compute tightened constraints
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            // Compute the standard deviation of the constraints
            terminalConstraintStdDev_(i) = std::sqrt(pointTransformation_->variance(terminalConstraintMatrix_.row(i)));

            // Tightened constraint
            out[i] = pointTransformation_->mean1D(terminalConstraintMatrix_.row(i)) + constraintTighteningCoeff_(numConstraints_ + i) * terminalConstraintStdDev_(i);
        }
    }

    void ResamplingGPProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        // Mapping of the output covariance
        Eigen::Map<Matrix> outCov(out.data() + numStates_, pointDim_, numStates_);

        // Get points and Cholesky decomposition of the covariance matrix that are both computed in hfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();
        const Matrix& covCholStateAndParam_ = llt_.matrixLLT();

        // Compute the derivative of the tightened constraint with respect to the constraints of the sigma points
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            terminalConstraintVec_.row(i) = (dmean1D_dPoints_ + pointTransformation_->dvar_dpoints(terminalConstraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(numConstraints_ + i) / 2.0 / std::max(stdDevMin_, terminalConstraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhTdx_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_, numStates_), 
                t, pointsInitial_.col(i), 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                terminalConstraintVec_.col(i),
                param
            );
            problemDescription_->dhTdp_vec(
                temp_vec_pointDim_numPoints_.segment(i*pointDim_+numStates_, numParams_), 
                t, 
                pointsInitial_.col(i), 
                pointsInitial_.col(i).segment(numStates_, numParams_), 
                terminalConstraintVec_.col(i),
                param
            );
        }
        out.segment(0, numStates_) = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);
        outCov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);
    }

    void ResamplingGPProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        out[0] = 0.0;

        // Get points that are computed in hTfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();

        // terminalConstraintVec_ is computed and saved in dhTdx_vec()
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint with respect to the prediction horizon
            problemDescription_->dhTdT_vec(tempScalar_, t, pointsInitial_.col(i), pointsInitial_.col(i).segment(numStates_, numParams_), terminalConstraintVec_.col(i), param);
            out[0] += tempScalar_(0);
        }
    }

    void ResamplingGPProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
    {
        const Matrix& covariance = state->covariance();

        // initial state
        x0_.segment(0, numStates_) = state->mean();
        for(typeInt i = 0; i < numStates_; ++i)
        {
            x0_.segment(numStates_ + i  * pointDim_, numStates_) = covariance.col(i);
            x0_.segment(2 * numStates_ + i * pointDim_, numParams_).setZero();
        }

        // initial parameter
        p0_ = param->mean();

        // set mean and covariance of the parameters
        meanStateAndParam_.bottomRows(numParams_) = param->mean();
        covStateAndParam_.block(numStates_, numStates_, numParams_, numParams_) = param->covariance();
    }

    void ResamplingGPProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        const Matrix& covariance = state->covariance();

        // initial state
        x0_.segment(0, numStates_) = state->mean();
        for(typeInt i = 0; i < numStates_; ++i)
        {
            x0_.segment(numStates_ + i  * pointDim_, numStates_) = covariance.col(i);
            x0_.segment(2 * numStates_ + i * pointDim_, numParams_).setZero();
        }
    }

    ctypeRNum* ResamplingGPProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* ResamplingGPProblemDescription::p0()
    {
        return p0_.data();
    }

    ResamplingGPProblemDescriptionPtr ResamplingGPProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                            PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec,
                                            const std::vector<typeInt>& dynamicsIndicesWithGP)
    {
        return ResamplingGPProblemDescriptionPtr(new ResamplingGPProblemDescription(problemDescription, constraintApproximation, pointTransformation, gaussianProcessVec, dynamicsIndicesWithGP));
    }

    ResamplingGPProblemDescriptionPtr ResamplingGPProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec,
                                            const std::vector<typeInt>& dynamicsIndicesWithGP)
    {
        return ResamplingGPProblemDescriptionPtr(new ResamplingGPProblemDescription(problemDescription, pointTransformation, gaussianProcessVec, dynamicsIndicesWithGP));
    }
}