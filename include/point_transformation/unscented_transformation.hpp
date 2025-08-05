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


#ifndef UNSCENTED_TRANSFORMATION_GENERATOR_HPP
#define UNSCENTED_TRANSFORMATION_GENERATOR_HPP

#include "point_transformation.hpp"
#include "util/grampc_s_util.hpp"

namespace grampc
{
    // Uncented transformation
    class UnscentedTransformation : public PointTransformation
    {
    public:
        // Constructor considering only the specified variables as uncertain
        UnscentedTransformation(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain);

        // Constructor considering all variables as uncertain
        UnscentedTransformation(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa);

        // Get points that represent a distribution
        virtual const Matrix& points(DistributionConstPtr dist) override;

        // Get points that represent a distribution with specified mean and Cholesky decomposition of the covariance matrix
        virtual const Matrix& points(VectorConstRef mean, MatrixConstRef covCholesky) override;

        // Get previously generated points
        virtual const Matrix& points() override;

        // Compute the mean of the distribution
        virtual const Vector& mean(MatrixConstRef points) override;

         // Compute the mean of the one-dimensional distribution
        virtual typeRNum mean1D(RowVectorConstRef points) override;

        // Compute the cross-covariance matrix of the random vectors x and y 
        virtual const Matrix& covariance(MatrixConstRef pointsX, MatrixConstRef pointsY) override;

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) override;

        // Jacobian dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec
        virtual const Vector& dmean_dpoints_vec(VectorConstRef vec) override;

        // Jacobian d(cov(X, Y) + cov(Y,X))/dX multiplied by vector vec, i.e. (d(vector(cov(X, Y) + cov(Y, X)))/dX)^T*vec
        virtual const Vector& dcov_dpointsX_vec(MatrixConstRef pointsY, VectorConstRef vec) override;

        // Jacobian d(cov(X, Y) + cov(Y,X))/dY multiplied by vector vec, i.e. (d(vector(cov(X, Y) + cov(Y, X)))/dY)^T*vec
        virtual const Vector& dcov_dpointsY_vec(MatrixConstRef pointsX, VectorConstRef vec) override;

        // Jacobian dpoints/dmean multiplied by vector vec, i.e. (dpoints/dmean)^T*vec
        virtual const Vector& dpoints_dmean_vec(VectorConstRef vec) override;

        // Jacobian dpoints/dcov multiplied by vector vec, i.e. (dpoints/vector(cov))^T*vec, covCholesky is the Cholesky decomposition of the covariance matrix
        virtual const Matrix& dpoints_dcov_vec(MatrixConstRef covCholesky, VectorConstRef vec) override;

        // dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() override;

        // Jacobian dvariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) override;

        // Return the number of points
        virtual typeInt numberOfPoints() const override;

    private:
        typeRNum alpha_;
        typeRNum beta_;
        typeRNum kappa_;
        typeInt dimX_;
        typeInt dimY_;
        typeInt numUncertainVariables_;
        typeInt numPoints_;
        Matrix normalizedPoints_;
        Matrix points_;
        Vector meanX_;
        Vector meanY_;
        Matrix covariance_;
        Matrix zeroMat_;
        Vector dmean_dpoints_vec_;
        Vector dcov_dpointsX_vec_;
        Vector dcov_dpointsY_vec_;
        Vector dvar_dpoints_;
        Vector dpoints_dmean_vec_;
        Matrix dpoints_dcov_vec_;
        Vector weightsMean_;
        Vector weightsVar_;
        Matrix diffToMeanX_;
        Matrix diffToMeanY_;
        RowVector diffToMean1D_;
        RowVector dCovY_;
        RowVector dCovX_;
        RowVector vecDiff_;
        Vector temp_vec_dimX_;
        Vector temp_vec_dimY_;
        Vector temp_vec_dimX_dimX_;
        typeRNum tempScalar;
    };

    // Constructor considering only the specified variables as uncertain
    PointTransformationPtr UT(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain);
    // Constructor considering all variables as uncertain
    PointTransformationPtr UT(typeInt dimX, typeInt dimY, typeRNum alpha, typeRNum beta, typeRNum kappa);
}

#endif // UNSCENTED_TRANSFORMATION_GENERATOR_HPP