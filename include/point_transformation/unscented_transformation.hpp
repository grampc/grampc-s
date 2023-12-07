#ifndef UNSCENTED_TRANSFORMATION_GENERATOR_HPP
#define UNSCENTED_TRANSFORMATION_GENERATOR_HPP

#include "point_transformation.hpp"

namespace grampc
{
    // Uncented transformation
    class UnscentedTransformation : public PointTransformation
    {
    public:
        // Constructor considering only the specified variables as uncertain
        UnscentedTransformation(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain);

        // Constructor considering all variables as uncertain
        UnscentedTransformation(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa);

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
        typeRNum alpha_;
        typeRNum beta_;
        typeRNum kappa_;
        typeInt dim_;
        std::vector<bool> considerUncertain_;
        typeInt numUncertainVariables_;
        typeInt numPoints_;
        Matrix normalizedPoints_;
        Matrix points_;
        Vector mean_;
        Matrix covariance_;
        Vector dmean_dpoints_vec_;
        Vector dcov_dpoints_vec_;
        Vector dvar_dpoints_;
        Vector weightsMean_;
        Vector weightsVar_;
        Matrix diffToMean_;
        RowVector diffToMean1D_;
        Matrix dcov_;
        Matrix diffToMean_jk_dx_;
        Matrix diffToMean_ik_dx_;
        typeRNum tempScalar;
    };

    PointTransformationPtr UT(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain);
    PointTransformationPtr UT(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa);
}

#endif // UNSCENTED_TRANSFORMATION_GENERATOR_HPP