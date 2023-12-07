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
        ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);

        // Constructor using a std vector
        ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder);

        // Constructor with identical quadrature order for all dimensions
        ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder);

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
                
        // All combinations of roots of the corresponding orthogonal polynomials
        const Matrix& roots();

        // All combinations of roots of the corresponding orthogonal polynomials
        const Vector& weights();

        // Normalized quadrature points for zero mean and covariance matrix = identity matrix
        const Matrix& normalizedPoints();

    private:
        typeInt dim_;
        typeInt numPoints_;
        Matrix normalizedPoints_;
        Matrix points_;
        Matrix roots_;
        Vector mean_;
        Matrix covariance_;
        Matrix covariance_temp_;
        Vector weights_;
        Vector dmean_dpoints_vec_;
        Vector dcov_dpoints_vec_;
        Vector dvar_dpoints_;
        Matrix diffToMean_;
        RowVector diffToMean1D_;
        Matrix dcov_;
        Matrix diffToMean_jk_dx_;
        Matrix diffToMean_ik_dx_;
        typeRNum tempScalar;
    };   

    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder);
    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder);
    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder);
}

#endif // COMPOSED_QUADRATURE_HPP