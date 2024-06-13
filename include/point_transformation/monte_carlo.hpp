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


#ifndef MONTE_CARLO_POINT_GENERATOR_HPP
#define MONTE_CARLO_POINT_GENERATOR_HPP

#include "point_transformation.hpp"

namespace grampc
{
    // Transformation using stochastic sampling points
    class MonteCarloTransformation : public PointTransformation
    {
    public:
        // Constructor for Monte-Carlo simulations with a fixed number of points
        MonteCarloTransformation(typeInt dimX, typeInt dimY, typeInt numberOfPoints, const RandomNumberGenerator& rng);

        // Get points that represent a distribution
        virtual const Matrix& points(DistributionConstPtr dist) override;

        // Get previously generated points
        virtual const Matrix& points() override;

        // Compute the mean of the distribution
        virtual const Vector& mean(MatrixConstRef points) override;

         // Compute the mean of the one-dimensional distribution
        virtual typeRNum mean1D(RowVectorConstRef points) override;

        // Compute the cross-covariance matrix of the random vectors x and y 
        virtual const Matrix& covariance(MatrixConstRef pointsX, MatrixConstRef pointsY) override;

        // Compute the covariance matrix of the random vector
        void covarianceMatrix(MatrixRef out, MatrixConstRef points);

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) override;

        // dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() override;

        // Jacobian dvariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) override;

        // Return the number of points
        virtual typeInt numberOfPoints() const override;

        // Set i-th point
        void setPoint(typeInt index, VectorConstRef newPoint);
    
    private:
        typeInt dimX_;
        typeInt dimY_;
        RandomNumberGenerator rng_;
        typeInt numberOfPoints_;
        Matrix points_;
        Vector meanX_;
        Vector meanY_;
        Matrix covariance_;
        Vector dmean1D_dpoints_;
        Vector dvar_dpoints_;
        typeRNum weight_;
        Matrix diffToMeanX_;
        Matrix diffToMeanY_;
        RowVector diffToMean1D_;
        typeRNum tempScalar;
    };

    // Constructor for Monte-Carlo simulations with a fixed number of points
    PointTransformationPtr MonteCarlo(typeInt dimX, typeInt dimY, typeInt numberOfPoints, const RandomNumberGenerator& rng);
}

#endif // MONTE_CARLO_POINT_GENERATOR_HPP