#define EIGEN_RUNTIME_NO_MALLOC
#include "problem_description/resampling_problem_description.hpp"

namespace grampc
{
    ResamplingProblemDescription::ResamplingProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                                               PointTransformationPtr pointTransformation)
        : ResamplingProblemDescription(problemDescription, constraintApproximation, pointTransformation, Matrix::Zero(1,1))
    {
        covWienerProcess_ = Matrix::Zero(numStates_, numStates_);
    }

    ResamplingProblemDescription::ResamplingProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                                               PointTransformationPtr pointTransformation, MatrixConstRef covWienerProcess)
        : numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient()),
          covWienerProcess_(covWienerProcess)
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
    }

    void ResamplingProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, Ng, &numConstraints_, NgT, &numTerminalConstraints_);

        *Nx = numStates_ + numStates_ * numStates_ + numStates_ * numParams_; // mean of states, covariance of states, cross-covariance of states and parameters
        *Nu = numControlInputs_;
        *Nh = numConstraints_;
        *Np = numParams_;
        *NhT = numTerminalConstraints_;
    }

    void ResamplingProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> cov(x + numStates_, pointDim_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Vector> d_stateMean(out, numStates_);
        Eigen::Map<Matrix> d_cov(out + numStates_, pointDim_, numStates_);
        
        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = stateMean;
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Transform points
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->ffct(pointsTransformed_.col(i).data(), t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_);
        }

        // Cross-covariance between initial points and transformed points
        const Matrix& crossCov = pointTransformation_->covariance(pointsInitial_, pointsTransformed_);

        // Time derivative of the mean
        d_stateMean = pointTransformation_->mean(pointsTransformed_);

        // Time derivative of the covariance
        d_cov = crossCov;
        d_cov.topRows(numStates_) += crossCov.topRows(numStates_).transpose() + covWienerProcess_;
    }

    void ResamplingProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {   
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> cov(x + numStates_, pointDim_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Vector> dd_stateMean(out, numStates_);
        Eigen::Map<Matrix> dd_cov(out + numStates_, pointDim_, numStates_);

         // Mapping of adj
        Eigen::Map<const Vector> adjMeanState(adj, numStates_);
        Eigen::Map<const Vector> adjCov(adj + numStates_, pointDim_ * numStates_);
        
        // covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = stateMean;
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Transform points
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->ffct(pointsTransformed_.col(i).data(), t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_);
        }

        // d(d_mean_state)/d(mean_state)
        const Vector& dmean_dpoints_vec = pointTransformation_->dmean_dpoints_vec(adjMeanState);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_, t, pointsInitial_.col(i).data(), dmean_dpoints_vec.data() + i * numStates_, u , pointsInitial_.col(i).data() + numStates_);
            problemDescription_->dfdp_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_ + numStates_, t, pointsInitial_.col(i).data(), dmean_dpoints_vec.data(), u , pointsInitial_.col(i).data() + numStates_);
        }
        dd_stateMean = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);

        // d(d_mean_state)/d(cov)
        dd_cov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);

        // d(d_cov)/d(mean)
        dcov_dPointsY_ = pointTransformation_->dcov_dpointsY_vec(pointsInitial_, adjCov);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_, t, pointsInitial_.col(i).data(), dcov_dPointsY_.data() + i * numStates_, u , pointsInitial_.col(i).data() + numStates_);
            problemDescription_->dfdp_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_ + numStates_, t, pointsInitial_.col(i).data(), dcov_dPointsY_.data(), u , pointsInitial_.col(i).data() + numStates_);
        }
        temp_vec_pointDim_numPoints_ += pointTransformation_->dcov_dpointsX_vec(pointsTransformed_, adjCov);
        dd_stateMean += pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);

        //d(d_cov)/d(cov)
        dd_cov += pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);

    }

    void ResamplingProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> cov(x + numStates_, pointDim_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Vector> outVec (out, numControlInputs_);
        outVec.setZero();

         // Mapping of adj
        Eigen::Map<const Vector> adjMeanState(adj, numStates_);
        Eigen::Map<const Vector> adjCov(adj + numStates_, pointDim_ * numStates_);
        
        // covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = stateMean;
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // d(d_mean)/du
        const Vector& dmean_dpoints_vec = pointTransformation_->dmean_dpoints_vec(adjMeanState);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(temp_vec_numInputs_.data(), t, pointsInitial_.col(i).data(), dmean_dpoints_vec.data() + i * numStates_, u , pointsInitial_.col(i).data() + numStates_);
            outVec += temp_vec_numInputs_;
        }

        // d(d_cov)/du
        dcov_dPointsY_ = pointTransformation_->dcov_dpointsY_vec(pointsInitial_, adjCov);
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(temp_vec_numInputs_.data(), t, pointsInitial_.col(i).data(), dcov_dPointsY_.data() + i * numStates_, u , pointsInitial_.col(i).data() + numStates_);
            outVec += temp_vec_numInputs_;
        }
    }

    void ResamplingProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->lfct(out, t, x, u, p, xdes, udes);
    }

    void ResamplingProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->dldx(out, t, x, u, p, xdes, udes);
    }

    void ResamplingProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->dldu(out, t, x, u, p, xdes, udes);
    }

    void ResamplingProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        problemDescription_->Vfct(out, t, x, p, xdes);
    }

    void ResamplingProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        problemDescription_->dVdx(out, t, x, p, xdes);
    }

    void ResamplingProblemDescription::dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
       problemDescription_->dVdT(out, t, x, p, xdes);
    }

    void ResamplingProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> cov(x + numStates_, pointDim_, numStates_);

        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = stateMean;
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(constraintMatrix_.data() + i * numConstraints_, t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_);
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

    void ResamplingProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the outputs
        Eigen::Map<Vector> outMean(out, numStates_);
        Eigen::Map<Matrix> outCov(out + numStates_, pointDim_, numStates_);

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
            problemDescription_->dhdx_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_, t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_, constraintVec_.col(i).data());
            problemDescription_->dhdp_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_ + numStates_, t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_, constraintVec_.col(i).data());
        }
        outMean = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);
        outCov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);
    }

    void ResamplingProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the output
        Eigen::Map<Vector> outVec(out, numControlInputs_);
        outVec.setZero();

        // Get points that are computed in hfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();

        // constraintVec is computed and saved in dhdx_vec()
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint with respect to the input
            problemDescription_->dhdu_vec(temp_vec_numInputs_.data(), t, pointsInitial_.col(i).data(), u, pointsInitial_.col(i).data() + numStates_, constraintVec_.col(i).data());
            outVec += temp_vec_numInputs_;
        }
    }

    void ResamplingProblemDescription::hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> cov(x + numStates_, pointDim_, numStates_);

        // Covariance matrix of states and parameter, bottom right corner is set in compute_x0_and_p0() and top right corner is not read
        covStateAndParam_.leftCols(numStates_) = cov;

        // Mean and Cholesky decomposition of the covariance matrix
        meanStateAndParam_.topRows(numStates_) = stateMean;
        const Matrix& covCholStateAndParam_ = llt_.compute(covStateAndParam_).matrixLLT();

        // Generate points
        const Matrix& pointsInitial_ = pointTransformation_->points(meanStateAndParam_, covCholStateAndParam_);

        // Evaluate constraints for each sigma point
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(terminalConstraintMatrix_.data() + i * numTerminalConstraints_, t, pointsInitial_.col(i).data(), pointsInitial_.col(i).data() + numStates_);
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

    void ResamplingProblemDescription::dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the outputs
        Eigen::Map<Vector> outMean(out, numStates_);
        Eigen::Map<Matrix> outCov(out + numStates_, pointDim_, numStates_);

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
            problemDescription_->dhTdx_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_, t, pointsInitial_.col(i).data(), pointsInitial_.col(i).data() + numStates_, terminalConstraintVec_.col(i).data());
            problemDescription_->dhTdp_vec(temp_vec_pointDim_numPoints_.data() + i * pointDim_ + numStates_, t, pointsInitial_.col(i).data(), pointsInitial_.col(i).data() + numStates_, terminalConstraintVec_.col(i).data());
        }
        outMean = pointTransformation_->dpoints_dmean_vec(temp_vec_pointDim_numPoints_);
        outCov = pointTransformation_->dpoints_dcov_vec(covCholStateAndParam_, temp_vec_pointDim_numPoints_);
    }

    void ResamplingProblemDescription::dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        out[0] = 0.0;

        // Get points that are computed in hTfct()
        const Matrix& pointsInitial_ = pointTransformation_->points();

        // terminalConstraintVec_ is computed and saved in dhTdx_vec()
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            // Derivative of the constraint with respect to the prediction horizon
            problemDescription_->dhTdT_vec(&tempScalar_, t, pointsInitial_.col(i).data(), pointsInitial_.col(i).data() + numStates_, terminalConstraintVec_.col(i).data());
            out[0] += tempScalar_;
        }
    }

    void ResamplingProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
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

    void ResamplingProblemDescription::compute_x0_and_p0(DistributionPtr state)
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

    ctypeRNum* ResamplingProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* ResamplingProblemDescription::p0()
    {
        return p0_.data();
    }
}