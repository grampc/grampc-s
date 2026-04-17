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


#include "problem_description/monte_carlo_problem_description.hpp"

namespace grampc
{
    MonteCarloProblemDescription::MonteCarloProblemDescription(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation)
        : ProblemDescription(
            problemDescription->Nx * pointTransformation->numberOfPoints(), 
            problemDescription->Nu, 
            problemDescription->Np * pointTransformation->numberOfPoints(),
            0,
            problemDescription->Nh * pointTransformation->numberOfPoints(),
            0,
            problemDescription->NhT * pointTransformation->numberOfPoints()
        ),
          numSigmaPoints_(pointTransformation->numberOfPoints()),
          problemDescription_(problemDescription),
          pointTransformation_(pointTransformation),
          tempScalar_(1)
    {
        // get the number of states, controls, parameters, and constraints of base problem description
        numStates_ = problemDescription_->Nx;
        numControlInputs_ = problemDescription_->Nu;
        numParams_ = problemDescription_->Np;
        numConstraints_ = problemDescription_->Nh;
        numTerminalConstraints_ = problemDescription_->NhT;

        // Distribution of states and parameters
        stateAndParam_ = MultiDist({Dist(numStates_), Dist(numParams_)});

        // allocate memory
        x0_ = Vector::Zero(numStates_ * numSigmaPoints_);
        p0_ = Vector::Zero(numParams_ * numSigmaPoints_);
        temp_vec_numInputs_ = Vector::Zero(numControlInputs_);
        temp_vec_numSigmaPoints_ = RowVector::Zero(numSigmaPoints_);
        dmean_dPoints_ = pointTransformation_->dmean1D_dpoints().transpose();
    }

    void MonteCarloProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
    {
        // f of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
           problemDescription_->ffct(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
        }
        // NOTE: diffusion of the samples (Wiener process noise) is only possible for discrete-time systems
        // or with a dedicated integrator (Euler-Maruyama integration) since the noise is not differentiable
    }

    void MonteCarloProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const GrampcParam& param)
    {
        // dfdx of all sigma points
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), adj.segment(i*numStates_, numStates_), param);
        }
    }

    void MonteCarloProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const GrampcParam& param)
    {
        // dfdu_vec of the first sigma point
        problemDescription_->dfdu_vec(out, t, x, adj, u, p, param);

        // Add dfdu_vec of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dfdu_vec(temp_vec_numInputs_, t, x.segment(i * numStates_, numStates_), u, p.segment(i * numParams_, numParams_), adj.segment(i * numStates_, numStates_), param);
            out += temp_vec_numInputs_;
        }
    }

    void MonteCarloProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->lfct(temp_vec_numSigmaPoints_.col(i), t, x.segment(numStates_, numStates_), u, p.segment(numParams_, numParams_), param);
        }

        // Output is the mean of the cost distribution
        out[0] = pointTransformation_->mean1D(temp_vec_numSigmaPoints_); 
    }

    void MonteCarloProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
    {
        // dl_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dldx(out.segment(j*numStates_, numStates_), t, x.segment(j*numStates_, numStates_), u, p.segment(j*numParams_, numParams_), param);
             out.segment(j*numStates_, numStates_) *= dmean_dPoints_(j);
        }
    }

    void MonteCarloProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
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

    void MonteCarloProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
    {
        // Cost of each sigmapoint
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->Vfct(temp_vec_numSigmaPoints_.col(i), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
        }

        // Output is the mean of the cost distribution
        out[0] = pointTransformation_->mean1D(temp_vec_numSigmaPoints_);
    }

    void MonteCarloProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
    {
        // dV_dx_ij = dmean_dpoints__j * dpoint_j/dx_ij
        for(typeInt j = 0; j < numSigmaPoints_; ++j)
        {
             problemDescription_->dVdx(out.segment(j*numStates_, numStates_), t, x.segment(j*numStates_, numStates_), p.segment(numParams_*j, numParams_), param);
             out.segment(j*numStates_, numStates_) *= dmean_dPoints_(j);
        }
    }

    void MonteCarloProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
    {
        // dVdT of the first sigma point
        problemDescription_->dVdT(out, t, x, p, param);
        out(0) *= dmean_dPoints_(0);
        
        // Add dVdT of the remaining sigma points
        for (int i = 1; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dVdT(tempScalar_, t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
            out(0) += tempScalar_(0) * dmean_dPoints_(i);
        }
    }

    void MonteCarloProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
    {
        // Inequality constraint of all samples
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hfct(out.segment(i*numConstraints_, numConstraints_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), param);
        }
    }

    void MonteCarloProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
    {   
        // dhdx of all samples
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), vec.segment(i*numConstraints_, numConstraints_), param);
        }
    }

    void MonteCarloProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
    {
        out.setZero();
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            /// Derivative of the constraint
            problemDescription_->dhdu_vec(temp_vec_numInputs_, t, x.segment(i*numStates_, numStates_), u, p.segment(i*numParams_, numParams_), vec.segment(i*numConstraints_, numConstraints_), param);
            out += temp_vec_numInputs_;
        }
    }

    void MonteCarloProblemDescription::hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
    {
        // Terminal inequality constraint of all samples
        for (typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->hTfct(out.segment(i*numTerminalConstraints_, numTerminalConstraints_), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), param);
        }
    }

    void MonteCarloProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
    {
        // dhTdx of all samples
        for (int i = 0; i < numSigmaPoints_; ++i)
        {
            problemDescription_->dhTdx_vec(out.segment(i*numStates_, numStates_), t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), vec.segment(i*numTerminalConstraints_, numTerminalConstraints_), param);
        }
    }

    void MonteCarloProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
    {
        out[0] = 0.0;
        for(typeInt i = 0; i < numSigmaPoints_; ++i)
        {
            /// Derivative of the constraint
            problemDescription_->dhTdT_vec(tempScalar_, t, x.segment(i*numStates_, numStates_), p.segment(i*numParams_, numParams_), vec.segment(numTerminalConstraints_*i, numTerminalConstraints_), param);
            out += tempScalar_;
        }
    }

    void MonteCarloProblemDescription::compute_x0_and_p0(DistributionPtr state, DistributionPtr param)
    {
        // Set combined distribution
        stateAndParam_->replaceDistribution(0, state);
        stateAndParam_->replaceDistribution(1, param);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_).reshaped(numStates_ * numSigmaPoints_, 1);
        p0_ = points.bottomRows(numParams_).reshaped(numParams_ * numSigmaPoints_, 1);
        // x0_ = points.topRows(numStates_);
        // p0_ = points.bottomRows(numParams_);
    }

    void MonteCarloProblemDescription::compute_x0_and_p0(DistributionPtr state)
    {
        // Set state distribution
        stateAndParam_->replaceDistribution(0, state);

        // Compute state and parameter points and copy them to matrices of the correct size
        const Matrix& points = pointTransformation_->points(stateAndParam_);
        x0_ = points.topRows(numStates_).reshaped(numStates_ * numSigmaPoints_, 1);
        p0_ = points.bottomRows(numParams_).reshaped(numParams_ * numSigmaPoints_, 1);
        // x0_ = points.topRows(numStates_);
        // p0_ = points.bottomRows(numParams_);
    }

    const Vector MonteCarloProblemDescription::x0()
    {
        return x0_;
    }

    const Vector MonteCarloProblemDescription::p0()
    {
        return p0_;
    }

    MonteCarloProblemDescriptionPtr MonteCarloProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation)
    {
        return MonteCarloProblemDescriptionPtr(new MonteCarloProblemDescription(problemDescription, pointTransformation));
    }
}
