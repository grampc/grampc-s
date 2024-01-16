#include "problem_description/taylor_problem_description.hpp"


namespace grampc
{
    TaylorProblemDescription::TaylorProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation)
        : TaylorProblemDescription(problemDescription, constraintApproximation, Matrix::Zero(1,1))
    {
        covWienerProcess_ = Matrix::Zero(numStates_, numStates_);
    }

    TaylorProblemDescription::TaylorProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef covWienerProcess)
        : problemDescription_(problemDescription),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient()),
          covWienerProcess_(covWienerProcess)
    {
        typeInt Ng, NgT;

        // Call ocp_dim to get the number of states, controls, parameters, and constraints
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, &Ng, &numConstraints_, &NgT, &numTerminalConstraints_);

        // allocate memory
        x0_ = Matrix::Zero(numStates_, 1 + numStates_ + numParams_);
        p0_ = Matrix::Zero(numParams_, 1 + numParams_);
        dfdx_ = Matrix::Zero(numStates_, numStates_);
        dfdp_ = Matrix::Zero(numStates_, numParams_);
        dfdxdx_ = Matrix::Zero(numStates_, numStates_ * numStates_);
        dfdxdp_ = Matrix::Zero(numStates_, numParams_ * numStates_);
        dfdxdu_ = Matrix::Zero(numStates_, numStates_ * numControlInputs_);
        dfdpdu_ = Matrix::Zero(numStates_, numParams_ * numControlInputs_);
        dhdx_ = Matrix::Zero(numConstraints_, numStates_);
        dhdvar_ = Vector::Zero(numConstraints_);
        dhdxdx_ = Matrix::Zero(numStates_, numStates_  * numConstraints_);
        dhdxdu_ = Matrix::Zero(numStates_, numControlInputs_ * numConstraints_);
        constraintStdDev_ = Vector::Zero(numConstraints_);
        dhTdx_ = Matrix::Zero(numTerminalConstraints_, numStates_);
        dhTdvar_ = Vector::Zero(numTerminalConstraints_);
        dhTdxdx_ = Matrix::Zero(numStates_, numStates_  * numTerminalConstraints_);
        dhTdxdT_ = Matrix::Zero(numStates_, numTerminalConstraints_);
        terminalConstraintStdDev_ = Vector::Zero(numTerminalConstraints_);
        temp_vec_numStates_ = Vector::Zero(numStates_);
        temp_mat_numStates_ = Matrix::Zero(numStates_, numStates_);
        temp_mat_numStates_numParams_ = Matrix::Zero(numStates_, numParams_);
        temp_mat_numStates_numControlInputs_ = Matrix::Zero(numStates_, numControlInputs_);
    }

    void TaylorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        problemDescription_->ocp_dim(&numStates_, &numControlInputs_, &numParams_, Ng, &numConstraints_, NgT, &numTerminalConstraints_);

        *Nx = numStates_ + numStates_ * numStates_ + numStates_ * numParams_; // mean of states, covariance of states, cross-covariance of states and parameters
        *Nu = numControlInputs_;
        *Nh = numConstraints_;
        *Np = numParams_ + numParams_ * numParams_;
        *NhT = numTerminalConstraints_;
    }

    void TaylorProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> stateCov(x + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Vector> paramMean(p, numParams_);
        Eigen::Map<const Matrix> paramCov(p + numParams_, numParams_, numParams_);

        // Mapping of the outputs
        Eigen::Map<Vector> d_stateMean(out, numStates_);
        Eigen::Map<Matrix> d_stateCov(out + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> d_stateParamCov(out + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdx(dfdx_.data(), t, stateMean.data(), u, paramMean.data());
        problemDescription_->dfdp(dfdp_.data(), t, stateMean.data(), u, paramMean.data());

        // d/dt mean = f(mean)
        problemDescription_->ffct(d_stateMean.data(), t, stateMean.data(), u, paramMean.data());

        // d/dt cov_state = dfdx(mean) * cov_state + cov_state * dfdx(mean)^T + dfdp(mean) * cov_state_param^T + cov_state_param * dfdp(mean)^T + covariance of Wiener process
        d_stateCov.noalias() = dfdx_ * stateCov + stateCov * dfdx_.transpose() + dfdp_ * stateParamCov.transpose() + stateParamCov * dfdp_.transpose() + covWienerProcess_;

        // d/dt cov_state_param = dfdx(mean) * cov_state_param + dfdp(mean) * cov_param
        d_stateParamCov.noalias() = dfdx_ * stateParamCov +  dfdp_ * paramCov;
    }

    void TaylorProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> stateCov(x + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Vector> paramMean(p, numParams_);
        Eigen::Map<const Matrix> paramCov(p + numParams_, numParams_, numParams_);

        // Mapping of the outputs
        Eigen::Map<Vector> dd_stateMean(out, numStates_);
        Eigen::Map<Matrix> dd_stateCov(out + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> dd_stateParamCov(out + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Mapping of adj
        Eigen::Map<const Matrix> adjCovState(adj + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> adjCovStateParam(adj + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdx(dfdx_.data(), t, stateMean.data(), u, paramMean.data());
        problemDescription_->dfdp(dfdp_.data(), t, stateMean.data(), u, paramMean.data());
        problemDescription_->dfdxdx(dfdxdx_.data(), t, stateMean.data(), u, paramMean.data());
        problemDescription_->dfdxdp(dfdxdp_.data(), t, stateMean.data(), u, paramMean.data());

        // d(d_mean_state)/d(mean_state)
        problemDescription_->dfdx_vec(dd_stateMean.data(), t, stateMean.data(), adj, u, paramMean.data());

        // d(d_cov_state)/d(cov_state)
        dd_stateCov.noalias() = dfdx_.transpose() * adjCovState + adjCovState * dfdx_;

        // d(d_cov_state)/d(cov_param_state)
        dd_stateParamCov.noalias() = adjCovState.transpose() * dfdp_ + adjCovState * dfdp_;

        // d(d_cov_state)/d(mean_state)
        for(typeInt i = 0; i < numStates_; ++i)
        {
            Eigen::Map<Matrix, 0, Eigen::OuterStride<>> dfdxdp_temp(dfdxdp_.data() + i * numStates_, numStates_, numParams_, Eigen::OuterStride<>(numStates_ * numStates_));
            temp_mat_numStates_.noalias() = dfdxdx_.block(0, i*numStates_, numStates_, numStates_) * stateCov +
                                             stateCov * dfdxdx_.block(0, i*numStates_, numStates_, numStates_).transpose() +
                                             dfdxdp_temp * stateParamCov.transpose() + stateParamCov * dfdxdp_temp.transpose();
            temp_mat_numStates_ = temp_mat_numStates_.cwiseProduct(adjCovState);
            dd_stateMean(i) += temp_mat_numStates_.sum();
        }
        
        // d(d_cov_param_state)/d(cov_param_state)
        dd_stateParamCov.noalias() += dfdx_.transpose() * adjCovStateParam;

        // d(d_cov_param_state)/d(mean_state)
        for(typeInt i = 0; i < numStates_; ++i)
        {
            Eigen::Map<Matrix, 0, Eigen::OuterStride<>> dfdxdp_temp(dfdxdp_.data() + i * numStates_, numStates_, numParams_, Eigen::OuterStride<>(numStates_ * numStates_));
            temp_mat_numStates_numParams_.noalias() = dfdxdx_.block(0, i*numStates_, numStates_, numStates_) * stateParamCov + dfdxdp_temp * paramCov;
            temp_mat_numStates_numParams_ = temp_mat_numStates_numParams_.cwiseProduct(adjCovStateParam);
            dd_stateMean(i) += temp_mat_numStates_numParams_.sum();
        }
    }

    void TaylorProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
    {
        // Mapping of the inputs
        Eigen::Map<const Vector> stateMean(x, numStates_);
        Eigen::Map<const Matrix> stateCov(x + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Vector> paramMean(p, numParams_);
        Eigen::Map<const Matrix> paramCov(p + numParams_, numParams_, numParams_);

        // Mapping of adj
        Eigen::Map<const Matrix> adjCovState(adj + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> adjCovStateParam(adj + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdxdu(dfdxdu_.data(), t, stateMean.data(), u, paramMean.data());
        problemDescription_->dfdpdu(dfdpdu_.data(), t, stateMean.data(), u, paramMean.data());

        // d(dmean_state)/du
        problemDescription_->dfdu_vec(out, t, stateMean.data(), adj, u, paramMean.data());

        // d(d_cov_state)/du
        for(typeInt i = 0; i < numControlInputs_; ++i)
        {
            temp_mat_numStates_.noalias() = dfdxdu_.block(0, i*numStates_, numStates_, numStates_) * stateCov +
                                            stateCov * dfdxdu_.block(0, i*numStates_, numStates_, numStates_).transpose() + 
                                            dfdpdu_.block(0, i*numParams_, numStates_, numParams_) * stateParamCov.transpose() +
                                            stateParamCov * dfdxdp_.block(0, i*numParams_, numStates_, numParams_).transpose();
            temp_mat_numStates_ = temp_mat_numStates_.cwiseProduct(adjCovState);
            out[i] += temp_mat_numStates_.sum();
        }

        // d(d_cov_param_state)/d(u)
        for(typeInt i = 0; i < numControlInputs_; ++i)
        {
            temp_mat_numStates_numParams_.noalias() = dfdxdu_.block(0, i*numStates_, numStates_, numStates_) * stateParamCov +
                                            dfdpdu_.block(0, i*numParams_, numStates_, numParams_) * paramCov;
            temp_mat_numStates_numParams_ = temp_mat_numStates_numParams_.cwiseProduct(adjCovStateParam);
            out[i] += temp_mat_numStates_numParams_.sum();
        }
    }

    void TaylorProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->lfct(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->dldx(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        problemDescription_->dldu(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        problemDescription_->Vfct(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        problemDescription_->dVdx(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        problemDescription_->dVdT(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // evaluate constraint at the mean
        problemDescription_->hfct(out, t, x, u, p);

        // evaluate constraint derivatives
        problemDescription_->dhdx(dhdx_.data(), t, x, u, p);

        // tighten constraints
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            temp_vec_numStates_.noalias() = covStates * dhdx_.row(i).transpose();

            constraintStdDev_(i) = std::sqrt(dhdx_.row(i) * temp_vec_numStates_);
            out[i] += constraintTighteningCoeff_(i) * constraintStdDev_(i);
        }
    }

    void TaylorProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Vector> outStateMean(out, numStates_);
        Eigen::Map<Matrix> outStateCov(out + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> outStateParamCov(out + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // derivative of the constraint
        problemDescription_->dhdxdx(dhdxdx_.data(), t, x, u, p);

        // Covariance between states and parameters does not effect the output
        outStateParamCov.setZero();

        // d(mean_h)/d(mean_x)
        problemDescription_->dhdx_vec(outStateMean.data(), t, x, u, p, vec);

        // d(var_h)/d(cov_x)
        outStateCov.setZero(); 
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            // use a lower bound of the standard deviation to prevent dividing by 0
            dhdvar_(i) = constraintTighteningCoeff_(i) * vec[i] / (2.0 * std::max(stdDevMin_, constraintStdDev_(i)));

            // dhdx_ is computed ind hfct
            outStateCov.noalias() += dhdvar_(i) * dhdx_.row(i).transpose() * dhdx_.row(i);
        }

        // d(var_h)/d(mean_x)
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            temp_mat_numStates_.noalias() = covStates * dhdxdx_.middleCols(i * numStates_, numStates_);

            // dhdx_ is computed in hfct
            outStateMean.noalias() += 2.0 * dhdvar_(i) * dhdx_.row(i) * temp_mat_numStates_;
        }
    }

    void TaylorProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // Mapping of the output
        Eigen::Map<Vector> outVec(out, numControlInputs_);

        // derivative of the constraint
        problemDescription_->dhdxdu(dhdxdu_.data(), t, x, u, p);

        // d(mean_h)/d(u)
        problemDescription_->dhdu_vec(outVec.data(), t, x, u, p, vec);

        // d(var_h)/d(u)
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            temp_mat_numStates_numControlInputs_.noalias() = covStates * dhdxdu_.middleCols(i * numControlInputs_, numControlInputs_);

            // dhdx_ is computed in hfct, dhdvar in dhdx_vec
            outVec.noalias() += 2.0 * dhdvar_(i) * dhdx_.row(i) * temp_mat_numStates_numControlInputs_;
        }
    }

    void TaylorProblemDescription::hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p)
    {
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // evaluate constraint at the mean
        problemDescription_->hTfct(out, t, x, p);

        // evaluate constraint derivatives
        problemDescription_->dhTdx(dhTdx_.data(), t, x, p);

        // tighten constraints
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            temp_vec_numStates_.noalias() = covStates * dhTdx_.row(i).transpose();

            terminalConstraintStdDev_(i) = std::sqrt(dhTdx_.row(i) * temp_vec_numStates_);
            out[i] += constraintTighteningCoeff_(numConstraints_ + i) * terminalConstraintStdDev_(i);
        }
    }

    void TaylorProblemDescription::dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Vector> outStateMean(out, numStates_);
        Eigen::Map<Matrix> outStateCov(out + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> outStateParamCov(out + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // derivative of the constraint
        problemDescription_->dhTdxdx(dhTdxdx_.data(), t, x, p);

        // Covariance between states and parameters does not effect the output
        outStateParamCov.setZero();

        // d(mean_hT)/d(mean_x)
        problemDescription_->dhTdx_vec(outStateMean.data(), t, x, p, vec);

        // d(var_hT)/d(cov_x)
        outStateCov.setZero(); 
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            // use a lower bound of the standard deviation to prevent dividing by 0
            dhTdvar_(i) = constraintTighteningCoeff_(numConstraints_ + i) * vec[i] / (2.0 * std::max(stdDevMin_, terminalConstraintStdDev_(i)));

            // dhTdx_ is computed ind hTfct
            outStateCov.noalias() += dhTdvar_(i) * dhTdx_.row(i).transpose() * dhTdx_.row(i);
        }

        // d(var_hT)/d(mean_x)
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            temp_mat_numStates_.noalias() = covStates * dhTdxdx_.middleCols(i * numStates_, numStates_);

            // dhTdx_ is computed in hTfct
            outStateMean.noalias() += 2.0 * dhTdvar_(i) * dhTdx_.row(i) * temp_mat_numStates_;
        }
    }

    void TaylorProblemDescription::dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec)
    {
        out[0] = 0.0;

        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x + numStates_, numStates_, numStates_);

        // d(mean_hT)/dT
        problemDescription_->dhTdT_vec(out, t, x, p, vec);

        // derivative of the constraint
        problemDescription_->dhTdxdT(dhTdxdT_.data(), t, x, p);

        // d(var_hT)/dT
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            temp_vec_numStates_.noalias() = covStates * dhTdxdT_.col(i);

            // dhTdx_ is computed in hfct, dhdvar in dhTdx_vec
            out[0] += 2.0 * dhTdvar_(i) * dhTdx_.row(i) * temp_vec_numStates_;
        }
    }

    void TaylorProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
    {
        x0_.col(0) = state->mean();
        x0_(Eigen::all, Eigen::seqN(1, numStates_)) = state->covariance();
        x0_(Eigen::all, Eigen::seqN(1+numStates_, numParams_)).setZero();

        p0_.col(0) = param->mean();
        p0_(Eigen::all, Eigen::seqN(1, numParams_)) = param->covariance();
    }

    void TaylorProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        x0_.col(0) = state->mean();
        x0_(Eigen::all, Eigen::seqN(1, numStates_)) = state->covariance();
        x0_(Eigen::all, Eigen::seqN(1+numStates_, numParams_)).setZero();
    }

    ctypeRNum* TaylorProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* TaylorProblemDescription::p0()
    {
        return p0_.data();
    }
}
