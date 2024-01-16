#include "problem_description/sigma_point_problem_description.hpp"

namespace grampc
{
    SigmaPointProblemDescription::SigmaPointProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                 PointTransformationPtr pointTransformation)
        : numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient())
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
    }

    void SigmaPointProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, Ng, &numConstraints_, NgT, &numTerminalConstraints_);
        *Nx = numSigmaPoints_ * numStates_;
        *Np = numSigmaPoints_ * numParams_;
        *Nu = numControlInputs_;
        *Nh = numConstraints_;
        *NhT = numTerminalConstraints_;
    }

    void SigmaPointProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // f of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
           problemDescription_->ffct(out + i * numStates_, t, x + i * numStates_, u, p + i * numParams_);
        }
    }

    void SigmaPointProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // dfdx of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(out + i * numStates_, t, x + i * numStates_, adj + i * numStates_, u, p + i * numParams_);
        }
    }

    void SigmaPointProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the output
        Eigen::Map<Vector> outVec(out, numControlInputs_);

        // dfdu_vec of the first sigma point
        problemDescription_->dfdu_vec(outVec.data(), t, x, adj, u, p);

        // Add dfdu_vec of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(temp_vec_numInputs_.data(), t, x + i * numStates_, adj + i * numStates_, u, p + i * numParams_);
            outVec += temp_vec_numInputs_;
        }
    }

    void SigmaPointProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->lfct(temp_vec_numSigmaPoints_.data() + i, t, x + i * numStates_, u, p + i * numParams_, xdes, udes);
        }

        // Output is the mean of the cost distribution
        *out = pointTransformation_->mean1D(temp_vec_numSigmaPoints_); 
    }

    void SigmaPointProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // Output Matrix
        Eigen::Map<Matrix> outMat(out, numStates_, numSigmaPoints_);

        // dl_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dldx(outMat.col(j).data(), t, x + j * numStates_, u, p + j * numParams_, xdes, udes);
             outMat.col(j) *= dmean_dPoints_(j);
        }
    }

    void SigmaPointProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // Output Vector
        Eigen::Map<Vector> outVec(out, numControlInputs_);

        // dldu of the first sigma point
        problemDescription_->dldu(temp_vec_numInputs_.data(), t, x, u, p, xdes, udes);
        outVec = temp_vec_numInputs_ * dmean_dPoints_(0);
        
        // Add dldu of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dldu(temp_vec_numInputs_.data(), t, x + i * numStates_, u, p + i * numParams_, xdes, udes);
            outVec += temp_vec_numInputs_ * dmean_dPoints_(i);
        }
    }

    void SigmaPointProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->Vfct(temp_vec_numSigmaPoints_.data() + i, t, x + i * numStates_, p + i * numParams_, xdes);
        }

        // Output is the mean of the cost distribution
        *out = pointTransformation_->mean1D(temp_vec_numSigmaPoints_);
    }

    void SigmaPointProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        // Output Matrix
        Eigen::Map<Matrix> outMat(out, numStates_, numSigmaPoints_);

        // dV_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dVdx(outMat.col(j).data(), t, x + j * numStates_, p + j * numParams_, xdes);
             outMat.col(j) *= dmean_dPoints_(j);
        }
    }

    void SigmaPointProblemDescription::dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        // dVdT of the first sigma point
        problemDescription_->dVdT(out, t, x, p, xdes);
        out[0] *= dmean_dPoints_(0);
        
        // Add dVdT of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dVdT(&tempScalar_, t, x + i * numStates_, p + i * numParams_, xdes);
            out[0] += tempScalar_ * dmean_dPoints_(i);
        }
    }

    void SigmaPointProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(constraintMatrix_.data() + i * numConstraints_, t, x + i * numStates_, u, p + i * numParams_);
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

    void SigmaPointProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {   
        // Compute the derivative of the tightened constraint with respect to the constraints of the sigma points
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            constraintVec_.row(i) = (dmean_dPoints_ + pointTransformation_->dvar_dpoints(constraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(i) / 2.0 / std::max(stdDevMin_, constraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhdx_vec(out + i * numStates_, t, x + i * numStates_, u, p + i * numParams_, constraintVec_.col(i).data());
        } 
    }

    void SigmaPointProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
       // Mapping of the output
        Eigen::Map<Vector> outMap(out, numControlInputs_);
        outMap.setZero();

        // constraintVec is computed and saved in dhdx_vec
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint
            problemDescription_->dhdu_vec(temp_vec_numInputs_.data(), t, x + i * numStates_, u, p + i * numParams_, constraintVec_.col(i).data());
            outMap += temp_vec_numInputs_;
        }
    }

    void SigmaPointProblemDescription::hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p)
    {
        // Evaluate the terminal constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(terminalConstraintMatrix_.data() + i * numTerminalConstraints_, t, x + i * numStates_, p + i * numParams_);
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

    void SigmaPointProblemDescription::dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        // Compute the derivative of the tightened constraint with respect to the terminal constraints of the sigma points
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            terminalConstraintVec_.row(i) = (dmean_dPoints_ + pointTransformation_->dvar_dpoints(terminalConstraintMatrix_.row(i)).transpose() * constraintTighteningCoeff_(numConstraints_ + i) / 2.0 / std::max(stdDevMin_, terminalConstraintStdDev_(i))) * vec[i];
        }

        // Derivative of the tightened constraints with respect to the states
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhTdx_vec(out + i * numStates_, t, x + i * numStates_, p + i * numParams_, terminalConstraintVec_.col(i).data());
        } 
    }

    void SigmaPointProblemDescription::dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        out[0] = 0.0;

        // terminalConstraintVec_ is computed and saved in dhTdx_vec
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint
            problemDescription_->dhTdT_vec(&tempScalar_, t, x + i * numStates_, p + i * numParams_, terminalConstraintVec_.col(i).data());
            out[0] += tempScalar_;
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

}
