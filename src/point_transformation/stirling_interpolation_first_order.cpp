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


#include "point_transformation/stirling_interpolation_first_order.hpp"

namespace grampc
{
    StirlingInterpolationFirstOrder::StirlingInterpolationFirstOrder(typeInt dimX, typeInt dimY, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    : dimX_(dimX),
      dimY_(dimY),
      stepSize_(stepSize),
      numUncertainVariables_(std::count(considerUncertain.begin(), considerUncertain.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dimX_, numPoints_),
      points_(dimX_, numPoints_),
      meanY_(dimY_),
      covariance_(dimX_, dimY_),
      dmean_dpoints_vec_(dimY_ * numPoints_),
      dmean1D_dpoints_vec_(numPoints_),
      dcov_dpointsX_vec_(Vector::Zero(dimX_ * numPoints_)),
      dcov_dpointsY_vec_(Vector::Zero(dimY_ * numPoints_)),
      dvar_dpoints_(numPoints_),
      dpoints_dmean_vec_(dimY_),
      dpoints_dcov_vec_(dimX_, dimY_),
      vecDiff_(dimX_ * dimX_),
      diffX_(dimX_, dimX_),
      diffY_(dimY, dimX),
      temp_vec_dimX_(dimX_),
      temp_vec_dimY_(dimY_),
      temp_vec_dimX_dimX_(dimX_ * dimX_)
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

        // Derivative of the one dimensional mean with respect to the points;
        dmean1D_dpoints_vec_ = Vector::Zero(numPoints_);
        dmean1D_dpoints_vec_(0) = 1.0;

        // Index for uncertain variables
        typeInt uncertIndex = 0;

        // Fill vector with indices of uncertain variables
        for(typeInt i = 0; i < dimX; ++i)
        {
            if(considerUncertain[i])
            {
                uncertainIndices[uncertIndex++] = i;
            }
        }

        // Set normalized Points
        Matrix A = Matrix::Zero(dimX, numUncertainVariables_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = stepSize_;
        }
        normalizedPoints_ << Vector::Zero(dimX), A, -A;
    }

    StirlingInterpolationFirstOrder::StirlingInterpolationFirstOrder(typeInt dimX, typeInt dimY, typeRNum stepSize)
    : StirlingInterpolationFirstOrder(dimX, dimY, stepSize, std::vector<bool>(dimX, true))
    {
    }

    const Matrix& StirlingInterpolationFirstOrder::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();

        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& StirlingInterpolationFirstOrder::points(VectorConstRef mean, MatrixConstRef covCholesky)
    {
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = mean;
            points_.col(i).noalias() += covCholesky * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& StirlingInterpolationFirstOrder::points()
    {
        return points_;
    }

    const Vector& StirlingInterpolationFirstOrder::mean(MatrixConstRef points)
    {
        meanY_ = points.col(0);
        return meanY_;
    }

    typeRNum StirlingInterpolationFirstOrder::mean1D(RowVectorConstRef points)
    {
        return points(0);
    }
    
    const Matrix& StirlingInterpolationFirstOrder::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        covariance_.setZero();
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            temp_vec_dimX_ = pointsX.col(i) - pointsX.col(i + numUncertainVariables_);
            temp_vec_dimY_ = pointsY.col(i) - pointsY.col(i + numUncertainVariables_);

            covariance_.noalias() += temp_vec_dimX_ * temp_vec_dimY_.transpose();
        }
        covariance_ /= (4.0 * stepSize_ * stepSize_);
        return covariance_;
    }

    typeRNum StirlingInterpolationFirstOrder::variance(RowVectorConstRef points)
    {
        tempScalar_ = 0.0;
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar2_ = points(i) - points(i + numUncertainVariables_);
            tempScalar_ += tempScalar2_ * tempScalar2_;
        }
        return tempScalar_ / (4.0 * stepSize_ * stepSize_);
    }

    const Vector& StirlingInterpolationFirstOrder::dmean_dpoints_vec(VectorConstRef vec)
    {
        for(typeInt i = 0; i < dimY_; ++i)
        {
            dmean_dpoints_vec_(i) = vec(i);
        }

        for(typeInt i = dimY_; i < dimY_*numPoints_; ++i)
        {
            dmean_dpoints_vec_(i) = 0.0;
        }

        return dmean_dpoints_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dmean1D_dpoints()
    {
        return dmean1D_dpoints_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dcov_dpointsX_vec(MatrixConstRef pointsY, VectorConstRef vec)
    {      
        // difference between points
        for(typeInt i = 1; i < dimX_ + 1; ++i)
        {
            diffY_.col(i-1) = pointsY.col(i) - pointsY.col(i + dimX_);
        }

        // Gradient of the first point is zero
        dcov_dpointsX_vec_.segment(0, dimX_).setZero();

        // Compute the gradient
        for(typeInt k = 1; k < dimX_ + 1; ++k)
        {
            for(typeInt i = 0; i < dimX_; ++i)
            {
                dcov_dpointsX_vec_(i + k * dimX_) = (diffY_.col(k-1).transpose() * vec(Eigen::seqN(i, dimY_, dimX_))).value() / (2.0 * stepSize_ * stepSize_);
                dcov_dpointsX_vec_(i + (k + dimX_) * dimX_) = - dcov_dpointsX_vec_(i + k * dimX_);
            }
        }
        return dcov_dpointsX_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dcov_dpointsY_vec(MatrixConstRef pointsX, VectorConstRef vec)
    {
        // difference between points
        for(typeInt i = 1; i < dimX_ + 1; ++i)
        {
            diffX_.col(i-1) = pointsX.col(i) - pointsX.col(i + dimX_);
        }

        // Gradient of the first point is zero
        dcov_dpointsY_vec_.segment(0, dimX_).setZero();

        // Compute the gradient
        for(typeInt k = 1; k < dimX_ + 1; ++k)
        {
            for(typeInt i = 0; i < dimY_; ++i)
            {
                dcov_dpointsY_vec_(i + k * dimY_) = (diffX_.col(k-1).transpose() * vec.segment(i * dimX_, dimX_)).value() / (2.0 * stepSize_ * stepSize_);
                dcov_dpointsY_vec_(i + (k + dimX_) * dimY_) = - dcov_dpointsY_vec_(i + k * dimY_);
            }
        }
        return dcov_dpointsY_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dpoints_dmean_vec(VectorConstRef vec)
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

    const Matrix& StirlingInterpolationFirstOrder::dpoints_dcov_vec(MatrixConstRef covCholesky, VectorConstRef vec)
    {
        // Vector of points 1 ... n minus Vector of points n+1 ... 2*n
        vecDiff_ = vec.segment(dimX_, dimX_ * dimX_).transpose() - vec.segment(dimX_ * (dimX_ + 1 ), dimX_ * dimX_).transpose();

        // upper half must be zero for the derivative of the Cholesky decompsition
        temp_vec_dimX_dimX_.setZero();
        for(typeInt k = 0; k < dimX_; ++k)
        {
            for(typeInt l = 0; l < dimY_; ++l)
            {
                // Compute derivative of the cholesky decomposition
                deriveCholesky(temp_vec_dimX_dimX_, covCholesky, k, l);

                // scale derivative 
                temp_vec_dimX_dimX_ *= stepSize_;

                // multiply derivative and vector
                dpoints_dcov_vec_(k, l) = vecDiff_ * temp_vec_dimX_dimX_;
            }
        }

        return dpoints_dcov_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dvar_dpoints(RowVectorConstRef points)
    {
        dvar_dpoints_(0) = 0.0;
        tempScalar_ = 1.0/(2.0 * stepSize_ * stepSize_);

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            dvar_dpoints_(i) = tempScalar_ * (points(i) - points(i + numUncertainVariables_));
            dvar_dpoints_(i + numUncertainVariables_) = -dvar_dpoints_(i);
        }
        return dvar_dpoints_;
    }

     typeInt StirlingInterpolationFirstOrder::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr StirlingFirstOrder(typeInt dimX, typeInt dimY, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new StirlingInterpolationFirstOrder(dimX, dimY, stepSize, considerUncertain));
    }

    PointTransformationPtr StirlingFirstOrder(typeInt dimX, typeInt dimY, typeRNum stepSize)
    {
        return PointTransformationPtr(new StirlingInterpolationFirstOrder(dimX, dimY, stepSize));
    }
}