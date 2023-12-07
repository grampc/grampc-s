#ifndef STIRLING_INTERPOLATION_SECOND_ORDER_HPP
#define STIRLING_INTERPOLATION_SECOND_ORDER_HPP

#include <iostream>
#include "point_transformation.hpp"

namespace grampc
{
    // Point transformaiton using Stirling's interpolation of second order
    class StirlingInterpolationSecondOrder : public PointTransformation
    {
    public:
        // Constructor considering only the specified variables as uncertain
        StirlingInterpolationSecondOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain);

        // Constructor considering all variables as uncertain
        StirlingInterpolationSecondOrder(typeInt dim, typeRNum stepSize);

        // Get points that represent a distribution
        virtual const Matrix& points(DistributionConstPtr dist) override;

        // Compute the mean of the distribution
        virtual const Vector& mean(MatrixConstRef points) override;

         // Compute the mean of the one-dimensional distribution
        virtual typeRNum mean1D(RowVectorConstRef points) override;

        // Compute the covariance matrix of the distribution
        virtual const Matrix& covariance(MatrixConstRef points) override;

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) override;

        // Jacobian dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec
        virtual const Vector& dmean_dpoints_vec(VectorConstRef vec) override;

        // dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() override;

        // Jacobian dcovariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec
        virtual const Vector& dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec) override;

        // Jacobian dvariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) override;

        // Return the number of points
        virtual typeInt numberOfPoints() const override;

    private:
        typeInt dim_;
        typeRNum stepsize_;
        typeRNum stepSizeSquared_;
        typeInt numUncertainVariables_;
        typeInt numPoints_;
        Matrix normalizedPoints_;
        Matrix points_;
        Vector weightsMean_;
        Vector weightsCov_;
        Vector mean_;
        Matrix covariance_;
        Vector dmean_dpoints_vec_;
        Vector dmean1D_dpoints_;
        Vector dcov_dpoints_vec_;
        Vector dvar_dpoints_;
        Vector tempVecDim_;
        Vector tempVecDim2_;
        typeRNum tempScalar_;
        typeRNum tempScalar2_;
        typeRNum tempScalar3_;
        typeRNum tempScalar4_;
    };   

    PointTransformationPtr StirlingSecondOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain);
    PointTransformationPtr StirlingSecondOrder(typeInt dim, typeRNum stepSize);
}

#endif // STIRLING_INTERPOLATION_SECOND_ORDER_HPP