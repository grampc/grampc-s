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


#include "point_transformation/stirling_interpolation_second_order.hpp"

namespace grampc
{
    StirlingInterpolationSecondOrder::StirlingInterpolationSecondOrder(typeInt dimX, typeInt dimY, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    : dimX_(dimX),
      dimY_(dimY),
      stepSize_(stepSize),
      stepSizeSquared_(stepSize*stepSize),
      numUncertainVariables_((typeInt) std::count(considerUncertain.begin(), considerUncertain.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dimX, numPoints_),
      points_(dimX, numPoints_),
      weightsMean_(2),
      weightsCov_(2),
      mean_(dimY),
      covariance_(dimX, dimY),
      dmean_dpoints_vec_(dimY*numPoints_),
      dcov_dpointsX_vec_(Vector::Zero(dimX_ * numPoints_)),
      dcov_dpointsY_vec_(Vector::Zero(dimY_ * numPoints_)),
      dvar_dpoints_(numPoints_),
      dpoints_dmean_vec_(dimY_),
      dpoints_dcov_vec_(dimX_, dimY_),
      dmean1D_dpoints_vec_(numPoints_),
      vecDiff_(dimX_ * dimX_),
      diffX_(dimX_, dimX_),
      diffY_(dimY, dimX),
      diffX2_(dimX_, dimX_),
      diffY2_(dimY, dimX),
      temp_vec_dimX_(dimX),
      temp_vec_dimY_(dimY),
      temp_vec_dimX_dimX_(dimX * dimX)
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

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

        // weights for mean and covariance computation
        weightsMean_(0) = (stepSizeSquared_ - numUncertainVariables_)/(stepSizeSquared_);
        weightsMean_(1) = 1.0 / (2.0 * stepSizeSquared_);
        weightsCov_(0) = 1.0 / (4.0 * stepSizeSquared_);
        weightsCov_(1) = (stepSizeSquared_ - 1.0) / (4.0 * stepSizeSquared_ * stepSizeSquared_);

        // Derivative of the one dimensional mean with respect to the points;
        dmean1D_dpoints_vec_(0) = weightsMean_(0);
        dmean1D_dpoints_vec_.segment(1, numPoints_ - 1).setConstant(weightsMean_(1));

        // Set normalized Points
        Matrix A = Matrix::Zero(dimX, numUncertainVariables_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = stepSize_;
        }
        normalizedPoints_ << Vector::Zero(dimX), A, -A;
    }

     StirlingInterpolationSecondOrder::StirlingInterpolationSecondOrder(typeInt dimX, typeInt dimY, typeRNum stepSize)
     : StirlingInterpolationSecondOrder(dimX, dimY, stepSize, std::vector<bool>(dimX, true))
    {
    }

    const Matrix& StirlingInterpolationSecondOrder::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }
        return points_;
    }

    const Matrix& StirlingInterpolationSecondOrder::points(VectorConstRef mean, MatrixConstRef covCholesky)
    {
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = mean;
            points_.col(i).noalias() += covCholesky * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& StirlingInterpolationSecondOrder::points()
    {
        return points_;
    }

    const Vector& StirlingInterpolationSecondOrder::mean(MatrixConstRef points)
    {
        mean_ = weightsMean_(0) * points.col(0);
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            mean_ += weightsMean_(1) * (points.col(i) + points.col(i + numUncertainVariables_));
        }
        return mean_;
    }

     typeRNum StirlingInterpolationSecondOrder::mean1D(RowVectorConstRef points)
    {
        tempScalar_ = weightsMean_(0) * points(0);
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar_ += weightsMean_(1) * (points(i) + points(i + numUncertainVariables_));
        }
        return tempScalar_;
    }

    const Matrix& StirlingInterpolationSecondOrder::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        covariance_.setZero();
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            temp_vec_dimX_ = pointsX.col(i) - pointsX.col(i + numUncertainVariables_);
            temp_vec_dimY_ = pointsY.col(i) - pointsY.col(i + numUncertainVariables_);

            covariance_.noalias() += weightsCov_(0) * temp_vec_dimX_ * temp_vec_dimY_.transpose();

            temp_vec_dimX_ = pointsX.col(i) + pointsX.col(i + numUncertainVariables_) - 2.0 * pointsX.col(0);
            temp_vec_dimY_ = pointsY.col(i) + pointsY.col(i + numUncertainVariables_) - 2.0 * pointsY.col(0);

            covariance_.noalias() += weightsCov_(1) * temp_vec_dimX_ * temp_vec_dimY_.transpose();
        }
        return covariance_;
    }

    typeRNum StirlingInterpolationSecondOrder::variance(RowVectorConstRef points)
    {
        tempScalar_ = 0.0;
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar2_ = points(i) - points(i + numUncertainVariables_);
            tempScalar3_ = points(i) + points(i + numUncertainVariables_) - 2.0 * points(0);
            tempScalar_ += tempScalar2_ * tempScalar2_ * weightsCov_(0) + tempScalar3_ * tempScalar3_ * weightsCov_(1);
        }
        return tempScalar_;
    }

     const Vector& StirlingInterpolationSecondOrder::dmean_dpoints_vec(VectorConstRef vec)
    {
        dmean_dpoints_vec_.segment(0, dimY_) = weightsMean_(0) * vec;

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            dmean_dpoints_vec_.segment(i*dimY_, dimY_) = weightsMean_(1) * vec;
            dmean_dpoints_vec_.segment((i + numUncertainVariables_) * dimY_, dimY_) = weightsMean_(1) * vec;
        }
        return dmean_dpoints_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dmean1D_dpoints()
    {
        return dmean1D_dpoints_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dcov_dpointsX_vec(MatrixConstRef pointsY, VectorConstRef vec)
    {      
        // difference between points
        for(typeInt i = 1; i < dimX_ + 1; ++i)
        {
            diffY_.col(i-1) = pointsY.col(i) - pointsY.col(i + dimX_);
            diffY2_.col(i-1) = pointsY.col(i) + pointsY.col(i + dimX_) - 2.0 * pointsY.col(0);
        }

        // Initialize gradient of the first point
        dcov_dpointsX_vec_.segment(0, dimX_).setZero();

        // Compute the gradient
        for(typeInt k = 1; k < dimX_ + 1; ++k)
        {
            for(typeInt i = 0; i < dimX_; ++i)
            {
                dcov_dpointsX_vec_(i + k * dimX_) = 2.0 * (diffY_.col(k-1).transpose() * vec(Eigen::seqN(i, dimY_, dimX_))).value() * weightsCov_(0);
                dcov_dpointsX_vec_(i + (k + dimX_) * dimX_) = - dcov_dpointsX_vec_(i + k * dimX_);
                                                  
                tempScalar_ = (diffY2_.col(k-1).transpose() * vec(Eigen::seqN(i, dimY_, dimX_))).value() * weightsCov_(1);
                dcov_dpointsX_vec_(i + k * dimX_) += tempScalar_;
                dcov_dpointsX_vec_(i + (k + dimX_) * dimX_) += tempScalar_;

                // Gradient of the first point
                dcov_dpointsX_vec_(i) += -2.0 * tempScalar_;
            }
        }
        return dcov_dpointsX_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dcov_dpointsY_vec(MatrixConstRef pointsX, VectorConstRef vec)
    {
        // difference between points
        for(typeInt i = 1; i < dimX_ + 1; ++i)
        {
            diffX_.col(i-1) = pointsX.col(i) - pointsX.col(i + dimX_);
            diffX2_.col(i-1) = pointsX.col(i) + pointsX.col(i + dimX_) - 2.0 * pointsX.col(0);
        }

        // Initialize gradient of the first point
        dcov_dpointsY_vec_.segment(0, dimY_).setZero();

        // Compute the gradient
        for(typeInt k = 1; k < dimX_ + 1; ++k)
        {
            for(typeInt i = 0; i < dimY_; ++i)
            {
                dcov_dpointsY_vec_(i + k * dimY_) =  2.0 * (diffX_.col(k-1).transpose() * vec.segment(i * dimX_, dimX_)).value() * weightsCov_(0);
                dcov_dpointsY_vec_(i + (k + dimX_) * dimY_) = - dcov_dpointsY_vec_(i + k * dimY_);
                                                  
                tempScalar_ = (diffX2_.col(k-1).transpose() * vec.segment(i * dimX_, dimX_)).value() * weightsCov_(1);
                dcov_dpointsY_vec_(i + k * dimY_) += tempScalar_;
                dcov_dpointsY_vec_(i + (k + dimX_) * dimY_) += tempScalar_;

                // Gradient of the first point
                dcov_dpointsY_vec_(i) += -2.0 * tempScalar_;
            }
        }
        return dcov_dpointsY_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dpoints_dmean_vec(VectorConstRef vec)
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

    const Matrix& StirlingInterpolationSecondOrder::dpoints_dcov_vec(MatrixConstRef covCholesky, VectorConstRef vec)
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

    const Vector& StirlingInterpolationSecondOrder::dvar_dpoints(RowVectorConstRef points)
    {
        tempScalar_ = 1.0/(2.0 * stepSize_ * stepSize_);
        tempScalar2_ = 2.0 * weightsCov_(1);

        dvar_dpoints_(0) = 0.0;

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar3_ = tempScalar2_ * (points(i) + points(i + numUncertainVariables_) - 2.0 * points(0));
            tempScalar4_ = tempScalar_ * (points(i) - points(i + numUncertainVariables_));
            dvar_dpoints_(i) = tempScalar3_ + tempScalar4_;
            dvar_dpoints_(i + numUncertainVariables_) = tempScalar3_ - tempScalar4_;
            dvar_dpoints_(0) -= 2.0 * tempScalar3_;
        }
        return dvar_dpoints_;
    }

     typeInt StirlingInterpolationSecondOrder::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr StirlingSecondOrder(typeInt dimX, typeInt dimY, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new StirlingInterpolationSecondOrder(dimX, dimY, stepSize, considerUncertain));
    }

    PointTransformationPtr StirlingSecondOrder(typeInt dimX, typeInt dimY, typeRNum stepSize)
    {
        return PointTransformationPtr(new StirlingInterpolationSecondOrder(dimX, dimY, stepSize));
    }
}