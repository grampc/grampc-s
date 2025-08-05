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


#ifndef POINT_TRANSFORMATION_HPP
#define POINT_TRANSFORMATION_HPP

#include "distribution/distribution.hpp"

namespace grampc
{
    // Transformation between a distribution and points that represent the distribution
    class PointTransformation
    {
    public:
        virtual ~PointTransformation() {}

        // Get points that represent a distribution
        virtual const Matrix& points(DistributionConstPtr dist) = 0;

        // Get points that represent a distribution with specified mean and Cholesky decomposition of the covariance matrix
        virtual const Matrix& points(VectorConstRef mean, MatrixConstRef covCholesky) {return zeroMatrix;}

        // Get previously generated points
        virtual const Matrix& points() = 0;

        // Compute the mean of the distribution
        virtual const Vector& mean(MatrixConstRef points) = 0;

         // Compute the mean of the one-dimensional distribution
        virtual typeRNum mean1D(RowVectorConstRef points) = 0;

        // Compute the cross-covariance matrix of the random vectors x and y 
        virtual const Matrix& covariance(MatrixConstRef pointsX, MatrixConstRef pointsY) = 0;

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) = 0;

        // Jacobian dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec
        virtual const Vector& dmean_dpoints_vec(VectorConstRef vec) {return zeroVector;}

        // Jacobian d(cov(X, Y) + cov(Y,X))/dX multiplied by vector vec, i.e. (d(vector(cov(X, Y) + cov(Y, X)))/dX)^T*vec
        virtual const Vector& dcov_dpointsX_vec(MatrixConstRef pointsY, VectorConstRef vec) {return zeroVector;}

        // Jacobian d(cov(X, Y) + cov(Y,X))/dY multiplied by vector vec, i.e. (d(vector(cov(X, Y) + cov(Y, X)))/dY)^T*vec
        virtual const Vector& dcov_dpointsY_vec(MatrixConstRef pointsX, VectorConstRef vec) {return zeroVector;}

        // Jacobian dpoints/dmean multiplied by vector vec, i.e. (dpoints/dmean)^T*vec
        virtual const Vector& dpoints_dmean_vec(VectorConstRef vec) {return zeroVector;}

        // Jacobian dpoints/dcov multiplied by vector vec, i.e. (dpoints/vector(cov))^T*vec, covCholesky is the Cholesky decomposition of the covariance matrix
        virtual const Matrix& dpoints_dcov_vec(MatrixConstRef covCholesky, VectorConstRef vec) {return zeroMatrix;}

        // dmean/dpoints for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() = 0;

        // Jacobian dvariance/dpoints for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) = 0;

        // Return the number of points
        virtual typeInt numberOfPoints() const = 0;

    private:
        Matrix zeroMatrix = Matrix::Zero(1, 1);
        Vector zeroVector = Vector::Zero(1);
    };

    // Alias
    typedef std::shared_ptr<PointTransformation> PointTransformationPtr;
    typedef std::shared_ptr<const PointTransformation> PointTransformationConstPtr;
}

#endif // POINT_TRANSFORMATION_HPP