#ifndef MONTE_CARLO_POINT_GENERATOR_HPP
#define MONTE_CARLO_POINT_GENERATOR_HPP

#include <algorithm>
#include "point_transformation.hpp"


namespace grampc
{
    // Transformation using stochastic sampling points
    class MonteCarloTransformation : public PointTransformation
    {
    public:
        // Constructor for Monte-Carlo simulations with a fixed number of points
        MonteCarloTransformation(typeInt dim, typeInt numberOfPoints, const RandomNumberGenerator& rng);

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
        RandomNumberGenerator rng_;
        typeInt numberOfPoints_;
        typeInt dim_;
        Matrix points_;
        Vector mean_;
        Matrix covariance_;
        Vector dmean_dpoints_vec_;
        Vector dmean1D_dpoints_;
        Vector dcov_dpoints_vec_;
        Vector dvar_dpoints_;
        typeRNum weight_;
        Matrix diffToMean_;
        RowVector diffToMean1D_;
        Matrix dcov_;
        Matrix diffToMean_jk_dx_;
        Matrix diffToMean_ik_dx_;
        typeRNum tempScalar;
    };

    PointTransformationPtr MonteCarlo(typeInt dim, typeInt numberOfPoints, const RandomNumberGenerator& rng);
}

#endif // MONTE_CARLO_POINT_GENERATOR_HPP