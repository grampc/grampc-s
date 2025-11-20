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


#include "point_transformation/unscented_transformation.hpp"

namespace grampc
{
    UnscentedTransformation::UnscentedTransformation(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain)
    : alpha_(alpha),
      beta_(beta),
      kappa_(kappa),
      dimX_(dimX),
      dimY_(dimY),
      numUncertainVariables_((typeInt) std::count(considerUncertain.begin(), considerUncertain.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dimX_, numPoints_),
      points_(dimX_, numPoints_),
      meanX_(dimX_),
      meanY_(dimY_),
      covariance_(dimX_, dimY_),
      zeroMat_(Matrix::Zero(dimX_, dimY_)),
      dmean_dpoints_vec_(dimY_ * numPoints_),
      dcov_dpointsX_vec_(Vector::Zero(dimX_ * numPoints_)),
      dcov_dpointsY_vec_(Vector::Zero(dimY_ * numPoints_)),
      dvar_dpoints_(numPoints_),
      dpoints_dmean_vec_(dimY_),
      dpoints_dcov_vec_(dimX_, dimY_),
      weightsMean_(numPoints_),
      weightsVar_(numPoints_),
      diffToMeanX_(dimX_, numPoints_),
      diffToMeanY_(dimY_, numPoints_),
      diffToMean1D_(numPoints_),
      dCovY_(dimY_),
      dCovX_(dimX_),
      vecDiff_(dimX_ * dimX_),
      temp_vec_dimX_(dimX_),
      temp_vec_dimY_(dimY_),
      temp_vec_dimX_dimX_(dimX_ * dimX_)
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

        // Index for uncertain variables
        typeInt uncertIndex = 0;

        // Fill vector with indices of uncertain variables
        for(typeInt i = 0; i < dimX_; ++i)
        {
            if(considerUncertain[i])
            {
                uncertainIndices[uncertIndex++] = i;
            }
        }

        // Set normalized Points
        Matrix A = Matrix::Zero(dimX_, numUncertainVariables_);
        tempScalar = alpha_ * std::sqrt(numUncertainVariables_ + kappa_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = tempScalar;
        }
        normalizedPoints_ << Vector::Zero(dimX_), A, -A;

        // Set weights
        ctypeRNum weight = 1.0 / (2.0 * alpha_*alpha_ * (numUncertainVariables_ + kappa_));
        weightsMean_ << 1.0 - numUncertainVariables_ / (alpha_*alpha_ * (numUncertainVariables_ + kappa_)), Vector::Constant(2*numUncertainVariables_, weight);
        weightsVar_ << weightsMean_(0) + 1.0 - alpha_*alpha_ + beta_, Vector::Constant(2*numUncertainVariables_, weight);
    }

    UnscentedTransformation::UnscentedTransformation(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa)
    : UnscentedTransformation(dimX, dimY, alpha, beta, kappa, std::vector<bool>(dimX, true))
    {
    }

    const Matrix& UnscentedTransformation::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();

        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& UnscentedTransformation::points(VectorConstRef mean, MatrixConstRef covCholesky)
    {
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = mean;
            points_.col(i).noalias() += covCholesky * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& UnscentedTransformation::points()
    {
        return points_;
    }

    const Vector& UnscentedTransformation::mean(MatrixConstRef points)
    {
        meanY_.noalias() = points * weightsMean_;
        return meanY_;
    }

    typeRNum UnscentedTransformation::mean1D(RowVectorConstRef points)
    {
        return points * weightsMean_;
    }

    const Matrix& UnscentedTransformation::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        // Compute mean of points
        meanX_.noalias() = pointsX * weightsMean_;
        meanY_.noalias() = pointsY * weightsMean_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMeanX_.col(i) = (pointsX.col(i) - meanX_) * weightsVar_(i);
            diffToMeanY_.col(i) = pointsY.col(i) - meanY_;
        }

        // Compute the covariance matrix
        covariance_.noalias() = diffToMeanX_ * diffToMeanY_.transpose();

        // improve numerical stability with this line
        //covariance_ = (covariance_.array().abs()<1e-12).select(zeroMat_,covariance_);
        return covariance_;
    }

    typeRNum UnscentedTransformation::variance(RowVectorConstRef points)
    {
        // mean
        tempScalar = points * weightsMean_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }

        // return variance
        return diffToMean1D_ * diffToMean1D_.transpose().cwiseProduct(weightsVar_);
    }

    const Vector& UnscentedTransformation::dmean_dpoints_vec(VectorConstRef vec)
    {
        // out = (w_1*v_1 w_1*v_2  w_1*v_3  ...  w_2*v_1  w_2*v_2   w_2*v_3  ...)^T

        // Compute and return output
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dmean_dpoints_vec_.segment(i * dimY_, dimY_) = weightsMean_(i) * vec;
        }

        return dmean_dpoints_vec_;
    }

    const Vector& UnscentedTransformation::dcov_dpointsX_vec(MatrixConstRef pointsY, VectorConstRef vec)
    {      
        // Compute mean of points
        meanY_.noalias() = pointsY * weightsMean_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMeanY_.col(i) = pointsY.col(i) - meanY_;
        }

        // weighted sum of all (Y_i - mean(Y)) --> required for the gradient below
        temp_vec_dimY_.noalias() = diffToMeanY_ * weightsVar_;

        // Compute the gradient:
        // dcov_ij/dx_kl = weightVar_k * y_kj - weightMean_k * sum_m(weightVar_m(y_mj - mean(y)_j))   if l = i
        // dcov_ij/dx_kl = 0 else
        for(typeInt k = 0; k < numPoints_; ++k)
        {
            dCovY_.noalias() = (weightsVar_(k) * diffToMeanY_.col(k) - weightsMean_(k) * temp_vec_dimY_).transpose();
            for(typeInt i = 0; i < dimX_; ++i)
            {
                // Multiplication with vec
                dcov_dpointsX_vec_(i + k * dimX_) = dCovY_ * (vec(Eigen::seqN(i, dimY_, dimX_)) + vec.segment(i * dimY_, dimY_));
            }
        }

        return dcov_dpointsX_vec_;
    }

    const Vector& UnscentedTransformation::dcov_dpointsY_vec(MatrixConstRef pointsX, VectorConstRef vec)
    {
        // Compute mean of points
        meanX_.noalias() = pointsX * weightsMean_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMeanX_.col(i) = pointsX.col(i) - meanX_;
        }

        // weighted sum of all (X_i - mean(X)) --> required for the gradient below
        temp_vec_dimX_.noalias() = diffToMeanX_ * weightsVar_;

        // Compute the gradient:
        // dcov_ij/dy_kl = weightVar_k * x_ki - weightMean_k * sum_m(weightVar_m(x_mi - mean(x)_i))   if l = j
        // dcov_ij/dy_kl = 0 else
        for(typeInt k = 0; k < numPoints_; ++k)
        {
            dCovX_.noalias() = (weightsVar_(k) * diffToMeanX_.col(k) - weightsMean_(k) * temp_vec_dimX_).transpose();
            for(typeInt i = 0; i < dimY_; ++i)
            {
                // Multiplication with vec
                dcov_dpointsY_vec_(i + k * dimY_) = dCovX_ * (vec.segment(i * dimX_, dimX_) + vec(Eigen::seqN(i, dimX_, dimY_)));
            }
        }

        return dcov_dpointsY_vec_;
    }

    const Vector& UnscentedTransformation::dpoints_dmean_vec(VectorConstRef vec)
    {
        for(typeInt i = 0; i < dimY_; ++i)
        {
            dpoints_dmean_vec_(i) = vec(i);
            for(typeInt j = 1; j < numPoints_; ++j)
            {
                dpoints_dmean_vec_(i) += vec(i + j * dimX_);
            }
        }

        return dpoints_dmean_vec_;
    }

    const Matrix& UnscentedTransformation::dpoints_dcov_vec(MatrixConstRef covCholesky, VectorConstRef vec)
    {
        // Vector of points 1 ... n minus Vector of points n+1 ... 2*n
        vecDiff_ = vec.segment(dimX_, dimX_ * dimX_).transpose() - vec.segment(dimX_ * (dimX_ + 1 ), dimX_ * dimX_).transpose();

        // scaling factor
        tempScalar = alpha_ * std::sqrt(dimX_ + kappa_);

        // upper half must be zero for the derivative of the Cholesky decompsition
        temp_vec_dimX_dimX_.setZero();
        for(typeInt k = 0; k < dimX_; ++k)
        {
            for(typeInt l = 0; l < dimY_; ++l)
            {
                // Compute derivative of the cholesky decomposition
                deriveCholesky(temp_vec_dimX_dimX_, covCholesky, k, l);

                // scale derivative 
                temp_vec_dimX_dimX_ *= tempScalar;

                // multiply derivative and vector
                dpoints_dcov_vec_(k, l) = vecDiff_ * temp_vec_dimX_dimX_;
            }
        }

        return dpoints_dcov_vec_;
    }

    const Vector& UnscentedTransformation::dmean1D_dpoints()
    {
        return weightsMean_;
    }

    const Vector& UnscentedTransformation::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points * weightsMean_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }

        tempScalar = - 2.0 * diffToMean1D_ * weightsMean_.cwiseProduct(weightsVar_); 
        
        // Derivative of the variance
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dvar_dpoints_(i) = 2.0 * weightsVar_(i) * diffToMean1D_(i) + tempScalar;
        }
        return dvar_dpoints_;
    }

    typeInt UnscentedTransformation::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr UT(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new UnscentedTransformation(dimX, dimY, alpha, beta, kappa, considerUncertain));
    }

    PointTransformationPtr UT(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa)
    {
        return PointTransformationPtr(new UnscentedTransformation(dimX, dimY, alpha, beta, kappa));
    }
    
}