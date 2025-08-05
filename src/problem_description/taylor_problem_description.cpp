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


#include "problem_description/taylor_problem_description.hpp"

namespace grampc
{
    TaylorProblemDescription::TaylorProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef diffMatrixWienerProcess)
        : problemDescription_(problemDescription),
          constraintTighteningCoeff_(constraintApproximation->tighteningCoefficient()),
          diffMatrixWienerProcess_(diffMatrixWienerProcess)
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

        // Array dimensions
        Nx_ = numStates_ + numStates_ * numStates_ + numStates_ * numParams_;
        Np_ = numParams_ + numParams_ * numParams_;
        Nu_ = numControlInputs_;
        Nh_ = numConstraints_;
        NhT_ = numTerminalConstraints_;
    }

    TaylorProblemDescription::TaylorProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation)
        : TaylorProblemDescription(problemDescription, constraintApproximation, Matrix::Zero(1,1))
    {
        diffMatrixWienerProcess_ = Matrix::Zero(numStates_, numStates_);
    }

    TaylorProblemDescription::TaylorProblemDescription(ProblemDescriptionPtr problemDescription, MatrixConstRef diffMatrixWienerProcess)
    : TaylorProblemDescription(problemDescription, Chebyshev(Vector(0)), diffMatrixWienerProcess)
    {
        if(numConstraints_ > 0)
        {
            std::cerr << "Constructor without chance constraint approximation is called but constraints exist in the problem description!" << std::endl;
        }
    }

    TaylorProblemDescription::TaylorProblemDescription(ProblemDescriptionPtr problemDescription)
    : TaylorProblemDescription(problemDescription, Chebyshev(Vector(0)), Matrix::Zero(1,1))
    {
        if(numConstraints_ > 0)
        {
            std::cerr << "Constructor without chance constraint approximation is called but constraints exist in the problem description!" << std::endl;
        }
        diffMatrixWienerProcess_ = Matrix::Zero(numStates_, numStates_);
    }

    void TaylorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx = Nx_;
        *Np = Np_;
        *Nu = Nu_;
        *Nh = Nh_;
        *NhT = NhT_;
        *Ng = 0;
        *NgT = 0;
    }

    void TaylorProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> stateCov(x.data() + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Matrix> paramCov(p.data() + numParams_, numParams_, numParams_);

        // Mapping of the outputs
        Eigen::Map<Matrix> d_stateCov(out.data() + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> d_stateParamCov(out.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdx(dfdx_.reshaped(), t, x, u, p);
        problemDescription_->dfdp(dfdp_.reshaped(), t, x, u, p);

        // d/dt mean = f(mean)
        problemDescription_->ffct(out, t, x, u, p);

        // d/dt cov_state = dfdx(mean) * cov_state + cov_state * dfdx(mean)^T + dfdp(mean) * cov_state_param^T + cov_state_param * dfdp(mean)^T + covariance of Wiener process
        d_stateCov.noalias() = dfdx_ * stateCov + stateCov * dfdx_.transpose() + dfdp_ * stateParamCov.transpose() + stateParamCov * dfdp_.transpose() + diffMatrixWienerProcess_;

        // d/dt cov_state_param = dfdx(mean) * cov_state_param + dfdp(mean) * cov_param
        d_stateParamCov.noalias() = dfdx_ * stateParamCov +  dfdp_ * paramCov;
    }

    void TaylorProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> stateCov(x.data() + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Matrix> paramCov(p.data() + numParams_, numParams_, numParams_);

        // Mapping of the outputs
        Eigen::Map<Matrix> dd_stateCov(out.data() + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> dd_stateParamCov(out.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Mapping of adj
        Eigen::Map<const Matrix> adjCovState(adj.data() + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> adjCovStateParam(adj.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdx(dfdx_.reshaped(), t, x, u, p);
        problemDescription_->dfdp(dfdp_.reshaped(), t, x, u, p);
        problemDescription_->dfdxdx(dfdxdx_.reshaped(), t, x, u, p);
        problemDescription_->dfdxdp(dfdxdp_.reshaped(), t, x, u, p);

        // d(d_mean_state)/d(mean_state)
        problemDescription_->dfdx_vec(out, t, x, adj, u, p);

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
            out(i) += temp_mat_numStates_.sum();
        }
        
        // d(d_cov_param_state)/d(cov_param_state)
        dd_stateParamCov.noalias() += dfdx_.transpose() * adjCovStateParam;

        // d(d_cov_param_state)/d(mean_state)
        for(typeInt i = 0; i < numStates_; ++i)
        {
            Eigen::Map<Matrix, 0, Eigen::OuterStride<>> dfdxdp_temp(dfdxdp_.data() + i * numStates_, numStates_, numParams_, Eigen::OuterStride<>(numStates_ * numStates_));
            temp_mat_numStates_numParams_.noalias() = dfdxdx_.block(0, i*numStates_, numStates_, numStates_) * stateParamCov + dfdxdp_temp * paramCov;
            temp_mat_numStates_numParams_ = temp_mat_numStates_numParams_.cwiseProduct(adjCovStateParam);
            out(i) += temp_mat_numStates_numParams_.sum();
        }
    }

    void TaylorProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> stateCov(x.data() + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> stateParamCov(x.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);
        Eigen::Map<const Matrix> paramCov(p.data() + numParams_, numParams_, numParams_);

        // Mapping of adj
        Eigen::Map<const Matrix> adjCovState(adj.data() + numStates_, numStates_, numStates_);
        Eigen::Map<const Matrix> adjCovStateParam(adj.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // Linearization of the system dynamics at the state mean and parameter mean
        problemDescription_->dfdxdu(dfdxdu_.reshaped(), t, x, u, p);
        problemDescription_->dfdpdu(dfdpdu_.reshaped(), t, x, u, p);

        // d(dmean_state)/du
        problemDescription_->dfdu_vec(out, t, x, adj, u, p);

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

    void TaylorProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        problemDescription_->lfct(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        problemDescription_->dldx(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
    {
        problemDescription_->dldu(out, t, x, u, p, xdes, udes);
    }

    void TaylorProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
    {
        problemDescription_->Vfct(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
    {
        problemDescription_->dVdx(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
    {
        problemDescription_->dVdT(out, t, x, p, xdes);
    }

    void TaylorProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
    {
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // evaluate constraint at the mean
        problemDescription_->hfct(out, t, x, u, p);

        // evaluate constraint derivatives
        problemDescription_->dhdx(dhdx_.reshaped(), t, x, u, p);

        // tighten constraints
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            temp_vec_numStates_.noalias() = covStates * dhdx_.row(i).transpose();

            constraintStdDev_(i) = std::sqrt(dhdx_.row(i) * temp_vec_numStates_);
            out[i] += constraintTighteningCoeff_(i) * constraintStdDev_(i);
        }
    }

    void TaylorProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Matrix> outStateCov(out.data() + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> outStateParamCov(out.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // derivative of the constraint
        problemDescription_->dhdxdx(dhdxdx_.reshaped(), t, x, u, p);

        // Covariance between states and parameters does not effect the output
        outStateParamCov.setZero();

        // d(mean_h)/d(mean_x)
        problemDescription_->dhdx_vec(out, t, x, u, p, vec);

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
            out.segment(0, numStates_).noalias() += 2.0 * dhdvar_(i) * dhdx_.row(i) * temp_mat_numStates_;
        }
    }

    void TaylorProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // derivative of the constraint
        problemDescription_->dhdxdu(dhdxdu_.reshaped(), t, x, u, p);

        // d(mean_h)/d(u)
        problemDescription_->dhdu_vec(out, t, x, u, p, vec);

        // d(var_h)/d(u)
        for(typeInt i = 0; i < numConstraints_; ++i)
        {
            temp_mat_numStates_numControlInputs_.noalias() = covStates * dhdxdu_.middleCols(i * numControlInputs_, numControlInputs_);

            // dhdx_ is computed in hfct, dhdvar in dhdx_vec
            out.noalias() += 2.0 * dhdvar_(i) * dhdx_.row(i) * temp_mat_numStates_numControlInputs_;
        }
    }

    void TaylorProblemDescription::hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p)
    {
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // evaluate constraint at the mean
        problemDescription_->hTfct(out, t, x, p);

        // evaluate constraint derivatives
        problemDescription_->dhTdx(dhTdx_.reshaped(), t, x, p);

        // tighten constraints
        for(typeInt i = 0; i < numTerminalConstraints_; ++i)
        {
            temp_vec_numStates_.noalias() = covStates * dhTdx_.row(i).transpose();

            terminalConstraintStdDev_(i) = std::sqrt(dhTdx_.row(i) * temp_vec_numStates_);
            out[i] += constraintTighteningCoeff_(numConstraints_ + i) * terminalConstraintStdDev_(i);
        }
    }

    void TaylorProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
    {
        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // Mapping of the outputs
        Eigen::Map<Matrix> outStateCov(out.data() + numStates_, numStates_, numStates_);
        Eigen::Map<Matrix> outStateParamCov(out.data() + numStates_ + numStates_ * numStates_, numStates_, numParams_);

        // derivative of the constraint
        problemDescription_->dhTdxdx(dhTdxdx_.reshaped(), t, x, p);

        // Covariance between states and parameters does not effect the output
        outStateParamCov.setZero();

        // d(mean_hT)/d(mean_x)
        problemDescription_->dhTdx_vec(out, t, x, p, vec);

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
            out.segment(0, numStates_).noalias() += 2.0 * dhTdvar_(i) * dhTdx_.row(i) * temp_mat_numStates_;
        }
    }

    void TaylorProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
    {
        out[0] = 0.0;

        // Mapping of the inputs
        Eigen::Map<const Matrix> covStates(x.data() + numStates_, numStates_, numStates_);

        // d(mean_hT)/dT
        problemDescription_->dhTdT_vec(out, t, x, p, vec);

        // derivative of the constraint
        problemDescription_->dhTdxdT(dhTdxdT_.reshaped(), t, x, p);

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
        x0_(Eigen::placeholders::all, Eigen::seqN(1, numStates_)) = state->covariance();
        x0_(Eigen::placeholders::all, Eigen::seqN(1+numStates_, numParams_)).setZero();

        p0_.col(0) = param->mean();
        p0_(Eigen::placeholders::all, Eigen::seqN(1, numParams_)) = param->covariance();
    }

    void TaylorProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        x0_.col(0) = state->mean();
        x0_(Eigen::placeholders::all, Eigen::seqN(1, numStates_)) = state->covariance();
        x0_(Eigen::placeholders::all, Eigen::seqN(1+numStates_, numParams_)).setZero();
    }

    ctypeRNum* TaylorProblemDescription::x0()
    {
        return x0_.data();
    }

    ctypeRNum* TaylorProblemDescription::p0()
    {
        return p0_.data();
    }

    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation)
    {
        return TaylorProblemDescriptionPtr(new TaylorProblemDescription(problemDescription, constraintApproximation));
    }

    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef diffMatrixWienerProcess)
    {
        return TaylorProblemDescriptionPtr(new TaylorProblemDescription(problemDescription, constraintApproximation, diffMatrixWienerProcess));
    }

    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription)
    {
        return TaylorProblemDescriptionPtr(new TaylorProblemDescription(problemDescription));
    }

    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, MatrixConstRef diffMatrixWienerProcess)
    {
        return TaylorProblemDescriptionPtr(new TaylorProblemDescription(problemDescription, diffMatrixWienerProcess));
    }
}
