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


#ifndef PCE_TRANSFORMATION_HPP
#define PCE_TRANSFORMATION_HPP

#include <map>
#include "composed_quadrature.hpp"
#include "polynomial/multivariate_polynomial.hpp"
#include "polynomial/orthogonal_polynomial_generator.hpp"

namespace grampc
{
    // Point transformation using polynomial chaos expansion and Gaussian quadrature
    class PCE_Transformation : public PointTransformation
    {
    public:
        // Constructor using an Eigen vector
        PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);
        
        // Constructor using a std vector
        PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const std::vector<typeInt>& quadratureOrder);

        // Constructor with identical quadrature order for all dimensions
        PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, typeInt quadratureOrder);

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

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) override;

        // dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() override;

        // Jacobian dvariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) override;

        // Return the number of points
        virtual typeInt numberOfPoints() const override;

    private:
        typeInt dimX_;
        typeInt dimY_;
        typeInt dimUncertain_;
        typeInt numPoints_;
        typeInt maxPolyOrder_;
        typeInt numPolynomials_;
        Matrix normalizedPoints_;
        Matrix points_;
        Vector PCE_coefficients_mean_;
        Matrix PCE_coefficientsX_;
        Matrix PCE_coefficientsY_;
        RowVector PCE_coefficients_var_;
        RowVector PCE_squared_coefficients_var_;
        Vector squaredNorms_;
        Vector squaredNorms_Var_;
        Vector temp_squaredNorms_Var_;
        Matrix covariance_;
        Vector weightsMean_;
        Matrix weightsCov_;
        Vector dvar_dpoints_vec_;
        Matrix dcoeff_dpoints;
        Matrix dcoeff_dpoints2;
    };

    // Constructor using an Eigen vector
    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);
    // Constructor using a std vector
    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const std::vector<typeInt>& quadratureOrder);
    // Constructor with identical quadrature order for all dimensions
    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, typeInt quadratureOrder);
}

#endif // PCE_TRANSFORMATION_HPP