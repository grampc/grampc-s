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


#ifndef COMPOSED_QUADRATURE_HPP
#define COMPOSED_QUADRATURE_HPP

#include "point_transformation.hpp"
#include "quadrature_rules/quadrature_rule.hpp"
#include "util/grampc_s_util.hpp"

namespace grampc
{
    // Quadrature rule for a multivariate distribution using the product rule for combination of the one-dimensional quadrature rules
    class ComposedQuadrature : public PointTransformation
    {
    public:
        // Constructor using an Eigen vector
        ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);

        // Constructor using a std vector
        ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder);

        // Constructor with identical quadrature order for all dimensions
        ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder);

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

        // dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() override;

        // Jacobian dvariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) override;

        // Return the number of points
        virtual typeInt numberOfPoints() const override;
                
        // All combinations of roots of the corresponding orthogonal polynomials
        const Matrix& roots();

        // All combinations of roots of the corresponding orthogonal polynomials
        const Vector& weights();

        // Normalized quadrature points for zero mean and covariance matrix = identity matrix
        const Matrix& normalizedPoints();

    private:
        typeInt dimX_;
        typeInt dimY_;
        typeInt numPoints_;
        Matrix normalizedPoints_;
        Matrix points_;
        Matrix roots_;
        Vector meanX_;
        Vector meanY_;
        Matrix covariance_;
        Vector weights_;
        Vector dvar_dpoints_;
        Matrix diffToMeanX_;
        Matrix diffToMeanY_;
        RowVector diffToMean1D_;
        typeRNum tempScalar;
    };   

    // Constructor using an Eigen vector
    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);
    // Constructor using a std vector
    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder);
    // Constructor with identical quadrature order for all dimensions
    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder);
}

#endif // COMPOSED_QUADRATURE_HPP