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


#include "problem_description/sigma_point_problem_description.hpp"

namespace grampc
{
    SigmaPointProblemDescription::SigmaPointProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                 PointTransformationPtr pointTransformation)
        : numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient()),
          tempScalar_(1)
    {
        typeInt Ng, NgT;

        // Call ocp_dim to get the number of states, controls, parameters, and constraints
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, &Ng, &numConstraints_, &NgT, &numTerminalConstraints_); 

        // Distribution of states and parameters
        stateAndParam_ = MultiDist({Dist(numStates_), Dist(numParams_)});

        // allocate memory
        x0_ = Matrix::Zero(numStates_, numSigmaPoints_);
        p0_ = Matrix::Zero(numParams_, numSigmaPoints_);
        temp_vec_numStates_ = Vector::Zero(numStates_);
        temp_vec_numInputs_ = Vector::Zero(numControlInputs_);
        temp_vec_numParams_ = Vector::Zero(numParams_);
        temp_vec_numSigmaPoints_ = RowVector::Zero(numSigmaPoints_);
        constraintMatrix_ = Matrix::Zero(numConstraints_, numSigmaPoints_);
        constraintStdDev_ = Vector::Zero(numConstraints_);
        constraintVec_ = Matrix::Zero(numConstraints_, numSigmaPoints_);
        terminalConstraintMatrix_ = Matrix::Zero(numTerminalConstraints_, numSigmaPoints_);
        terminalConstraintStdDev_ = Vector::Zero(numTerminalConstraints_);
        terminalConstraintVec_ = Matrix::Zero(numTerminalConstraints_, numSigmaPoints_);
        dmean_dPoints_ = pointTransformation_->dmean1D_dpoints().transpose();

        // Array dimensions
        Nx_ = numSigmaPoints_ * numStates_;
        Np_ = numSigmaPoints_ * numParams_;
        Nu_ = numControlInputs_;
        Nh_ = numConstraints_;
        NhT_ = numTerminalConstraints_; 
    }

    SigmaPointProblemDescription::SigmaPointProblemDescription(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation)
    : SigmaPointProblemDescription(problemDescription, Chebyshev(Vector(0)), pointTransformation)
    {
        if(numConstraints_ > 0)
        {
            std::cerr << "Constructor without chance constraint approximation is called but constraints exist in the problem description!" << std::endl;
        }
    }

    void SigmaPointProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Np = Np_;
        *Nu = Nu_;
        *Nh = Nh_;
        *NhT = NhT_;
        *Ng = 0;
        *NgT = 0;
    }

    void SigmaPointProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // f of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
           problemDescription_->ffct(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
        }
    }

    void SigmaPointProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param)
    {
        // dfdx of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), adj.segment(i*numStates_, numStates_), param);
        }
    }

    void SigmaPointProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param)
    {
        // dfdu_vec of the first sigma point
        problemDescription_->dfdu_vec(out, t, x, u, p, adj, param);

        // Add dfdu_vec of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(temp_vec_numInputs_, t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), adj.segment(i*numStates_, numStates_), param);
            out += temp_vec_numInputs_;
        }
    }

    void SigmaPointProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->lfct(temp_vec_numSigmaPoints_.col(i), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
        }

        // Output is the mean of the cost distribution
        out(0) = pointTransformation_->mean1D(temp_vec_numSigmaPoints_); 
    }

    void SigmaPointProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Output Matrix
        Eigen::Map<Matrix> outMat(out.data(), numStates_, numSigmaPoints_);

        // dl_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dldx(outMat.col(j), t, x.segment(j*numStates_, numStates_), u, p.segment(j*numParams_, numParams_), param);
             outMat.col(j) *= dmean_dPoints_(j);
        }
    }

    void SigmaPointProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // dldu of the first sigma point
        problemDescription_->dldu(temp_vec_numInputs_, t, x, u, p, param);
        out = temp_vec_numInputs_ * dmean_dPoints_(0);
        
        // Add dldu of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dldu(temp_vec_numInputs_, t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
            out += temp_vec_numInputs_ * dmean_dPoints_(i);
        }
    }

    void SigmaPointProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->Vfct(temp_vec_numSigmaPoints_.col(i), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
        }

        // Output is the mean of the cost distribution
        out(0) = pointTransformation_->mean1D(temp_vec_numSigmaPoints_);
    }

    void SigmaPointProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Output Matrix
        Eigen::Map<Matrix> outMat(out.data(), numStates_, numSigmaPoints_);

        // dV_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dVdx(outMat.col(j), t, x.segment(j*numStates_, numStates_), p.segment(j*numParams_, numParams_), param);
             outMat.col(j) *= dmean_dPoints_(j);
        }
    }

    void SigmaPointProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // dVdT of the first sigma point
        problemDescription_->dVdT(out, t, x, p, param);
        out[0] *= dmean_dPoints_(0);
        
        // Add dVdT of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dVdT(tempScalar_, t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
            out[0] += tempScalar_(0) * dmean_dPoints_(i);
        }
    }

    void SigmaPointProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(constraintMatrix_.col(i), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
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

    void SigmaPointProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {   
        // Compute the derivative of the tightened constraint with respect to the constraints of the sigma points
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            constraintVec_.row(i) = (dmean_dPoints_ + pointTransformation_->dvar_dpoints(constraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(i) / 2.0 / std::max(stdDevMin_, constraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), constraintVec_.col(i), param);
        } 
    }

    void SigmaPointProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        out.setZero();

        // constraintVec is computed and saved in dhdx_vec
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint
            problemDescription_->dhdu_vec(temp_vec_numInputs_, t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), constraintVec_.col(i), param);
            out += temp_vec_numInputs_;
        }
    }

    void SigmaPointProblemDescription::hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
    {
        // Evaluate the terminal constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(terminalConstraintMatrix_.col(i), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
        }
            
        // Compute tightened constraints
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            // Compute the standard deviation of the terminal constraints
            terminalConstraintStdDev_(i) = std::sqrt(pointTransformation_->variance(terminalConstraintMatrix_.row(i)));

            // Tightened constraint
            out[i] = pointTransformation_->mean1D(terminalConstraintMatrix_.row(i)) + constraintTighteningCoeff_(numConstraints_ + i) * terminalConstraintStdDev_(i);
        }
    }

    void SigmaPointProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        // Compute the derivative of the tightened constraint with respect to the terminal constraints of the sigma points
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            terminalConstraintVec_.row(i) = (dmean_dPoints_ + pointTransformation_->dvar_dpoints(terminalConstraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(numConstraints_ + i) / 2.0 / std::max(stdDevMin_, terminalConstraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhTdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), terminalConstraintVec_.col(i), param);
        } 
    }

    void SigmaPointProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
    {
        out[0] = 0.0;

        // terminalConstraintVec_ is computed and saved in dhTdx_vec
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint
            problemDescription_->dhTdT_vec(tempScalar_, t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), terminalConstraintVec_.col(i), param);
            out[0] += tempScalar_(0);
        }
    }

    void SigmaPointProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
    {
        // Set combined distribution
        stateAndParam_->replaceDistribution(0, state);
        stateAndParam_->replaceDistribution(1, param);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_);
        p0_ = points.bottomRows(numParams_);
    }

    void SigmaPointProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        // Set combined distribution
        stateAndParam_->replaceDistribution(0, state);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_);
        p0_ = points.bottomRows(numParams_);
    }

    ctypeRNum* SigmaPointProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* SigmaPointProblemDescription::p0()
    {
        return p0_.data();
    }

    SigmaPointProblemDescriptionPtr SigmaPointProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation)
    {
        return SigmaPointProblemDescriptionPtr(new SigmaPointProblemDescription(problemDescription, constraintApproximation, pointTransformation));
    }

    SigmaPointProblemDescriptionPtr SigmaPointProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation)
    {
        return SigmaPointProblemDescriptionPtr(new SigmaPointProblemDescription(problemDescription, pointTransformation));
    }
}
