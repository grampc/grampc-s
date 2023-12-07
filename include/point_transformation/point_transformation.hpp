#ifndef POINT_TRANSFORMATION_HPP
#define POINT_TRANSFORMATION_HPP

#include <vector>
#include <memory>
#include "distribution/distribution.hpp"

namespace grampc
{
    // Transformation between a distribution and points that represent the distribution
    class PointTransformation
    {
    public:
        virtual ~PointTransformation() {};

        // Get points that represent a distribution
        virtual const Matrix& points(DistributionConstPtr dist) = 0;

        // Compute the mean of the distribution
        virtual const Vector& mean(MatrixConstRef points) = 0;

         // Compute the mean of the one-dimensional distribution
        virtual typeRNum mean1D(RowVectorConstRef points) = 0;

        // Compute the covariance matrix of the distribution
        virtual const Matrix& covariance(MatrixConstRef points) = 0;

        // Compute the variance of the one-dimensional distribution
        virtual typeRNum variance(RowVectorConstRef points) = 0;

        // Jacobian dmean/dpoints multiplied by vector vec, i.e. (dmean/dpoints)^T*vec
        virtual const Vector& dmean_dpoints_vec(VectorConstRef vec) = 0;

        // dmean/dpoints for a one-dimensional distribution
        virtual const Vector& dmean1D_dpoints() = 0;

        // Jacobian dcovariance/dpoints multiplied by vector vec, i.e. (d(vector(covariance))/dpoints)^T*vec
        virtual const Vector& dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec) = 0;

        // Jacobian dvariance/dpoints for a one-dimensional distribution
        virtual const Vector& dvar_dpoints(RowVectorConstRef points) = 0;

        // Return the number of points
        virtual typeInt numberOfPoints() const = 0;
    };

    // Alias
    typedef std::shared_ptr<PointTransformation> PointTransformationPtr;
    typedef std::shared_ptr<const PointTransformation> PointTransformationConstPtr;
}

#endif // POINT_TRANSFORMATION_HPP