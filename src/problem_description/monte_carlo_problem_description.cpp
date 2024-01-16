#include "problem_description/monte_carlo_problem_description.hpp"

namespace grampc
{

    MonteCarloProblemDescription::MonteCarloProblemDescription(StochasticProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation)
        : numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation)
    {
        typeInt Ng, NgT;

        // Call ocp_dim to get the number of states, controls, parameters, and constraints
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, &Ng, &numConstraints_, &NgT, &numTerminalConstraints_);

        // Distribution of states and parameters
        stateAndParam_ = MultiDist({Dist(numStates_), Dist(numParams_)});

        // allocate memory
        x0_ = Matrix::Zero(numStates_, numSigmaPoints_);
        p0_ = Matrix::Zero(numParams_, numSigmaPoints_);
        temp_vec_numInputs_ = Vector::Zero(numControlInputs_);
        temp_vec_numSigmaPoints_ = RowVector::Zero(numSigmaPoints_);
        dmean_dPoints_ = pointTransformation_->dmean1D_dpoints().transpose();
    }

    void MonteCarloProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, Ng, &numConstraints_, NgT, &numTerminalConstraints_);
        *Nx = numSigmaPoints_ * numStates_;
        *Np = numSigmaPoints_ * numParams_;
        *Nu = numControlInputs_;
        *Nh = numConstraints_ * numSigmaPoints_;
        *NhT = numTerminalConstraints_ * numSigmaPoints_;
    }

    void MonteCarloProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // f of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
           problemDescription_->ffct(out + i * numStates_, t, x + i * numStates_, u, p + i * numParams_);
        }
        // NOTE: diffusion of the samples (Wiener process noise) is only possible for discrete-time systems
        // or with a dedicated integrator (Euler-Maruyama integration) since the noise is not differentiable
    }

    void MonteCarloProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // dfdx of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(out + i * numStates_, t, x + i * numStates_, adj + i * numStates_, u, p + i * numParams_);
        }
    }

    void MonteCarloProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
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

    void MonteCarloProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->lfct(temp_vec_numSigmaPoints_.data() + i, t, x + i * numStates_, u, p + i * numParams_, xdes, udes);
        }

        // Output is the mean of the cost distribution
        *out = pointTransformation_->mean1D(temp_vec_numSigmaPoints_); 
    }

    void MonteCarloProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
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

    void MonteCarloProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
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

    void MonteCarloProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->Vfct(temp_vec_numSigmaPoints_.data() + i, t, x + i * numStates_, p + i * numParams_, xdes);
        }

        // Output is the mean of the cost distribution
        *out = pointTransformation_->mean1D(temp_vec_numSigmaPoints_);
    }

    void MonteCarloProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
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

    void MonteCarloProblemDescription::dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
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

    void MonteCarloProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // Inequality constraint of all samples
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(out + i * numConstraints_, t, x + i * numStates_, u, p + i * numParams_);
        }
    }

    void MonteCarloProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {   
        // dhdx of all samples
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhdx_vec(out + i * numStates_, t, x + i * numStates_, u, p + i * numParams_, vec + i * numConstraints_);
        }
    }

    void MonteCarloProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
       // Mapping of the output
        Eigen::Map<Vector> outMap(out, numControlInputs_);
        outMap.setZero();

        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            /// Derivative of the constraint
            problemDescription_->dhdu_vec(temp_vec_numInputs_.data(), t, x + i * numStates_, u, p + i * numParams_, vec + i * numConstraints_);
            outMap += temp_vec_numInputs_;
        }
    }

    void MonteCarloProblemDescription::hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p)
    {
        // Terminal inequality constraint of all samples
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(out + i * numTerminalConstraints_, t, x + i * numStates_, p + i * numParams_);
        }
    }

    void MonteCarloProblemDescription::dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        // dhTdx of all samples
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhTdx_vec(out + i * numStates_, t, x + i * numStates_, p + i * numParams_, vec + i * numTerminalConstraints_);
        }
    }

    void MonteCarloProblemDescription::dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        out[0] = 0.0;
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            /// Derivative of the constraint
            problemDescription_->dhTdT_vec(&tempScalar_, t, x + i * numStates_, p + i * numParams_, vec + i * numTerminalConstraints_);
            out[0] += tempScalar_;
        }
    }

    void MonteCarloProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
    {
        // Set combined distribution
        stateAndParam_->replaceDistribution(0, state);
        stateAndParam_->replaceDistribution(1, param);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_);
        p0_ = points.bottomRows(numParams_);
    }

    void MonteCarloProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        // Set state distribution
        stateAndParam_->replaceDistribution(0, state);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_);
        p0_ = points.bottomRows(numParams_);
    }

    ctypeRNum* MonteCarloProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* MonteCarloProblemDescription::p0()
    {
        return p0_.data();
    }

}
