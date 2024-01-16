#ifndef DISTRIBUTION_HPP
#define DISTRIBUTION_HPP

#include <vector>
#include <memory>
#include <iostream>
#include "problem_description/grampc_interface.hpp"
#include "util/grampc_s_constants.hpp"

namespace grampc
{
    // Probability distribution with corresponding polynomial family for Gaussian quadrature
    class Distribution
    {
    public:
        virtual ~Distribution() {};

        // Constructor for a distribution of specified dimension and unspecified attributes 
        Distribution(typeInt dim);

        // Constructor for a distribution with mean, covariance matrix, and polynomial family for Gaussian quadrature
        Distribution(VectorConstRef mean, MatrixConstRef covariance, const std::vector<PolynomialFamily>& poly);

        // Constructor for a distribution with mean, covariance matrix, Cholesky decomposition of the covariance matrix, and polynomial family for Gaussian quadrature
        Distribution(VectorConstRef mean, MatrixConstRef covariance, MatrixConstRef covChol,  const std::vector<PolynomialFamily>& poly);

        // Get the mean vector of the distribution
        const Vector& mean() const;

        // Get the covariance matrix of the distribution
        const Matrix& covariance() const;

        // Get the Cholesky decomposition of the covariance matrix
        const Matrix& covCholesky() const;

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const;

        // Replace the mean and the covariance of this distribution
        void replaceMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance);

        // Get the diminsion of the distribution
        typeInt dimension() const;

        // Return vector of polynomial families
        const std::vector<PolynomialFamily>& polynomialFamily() const;

        // Set new mean vector and covariance matrix of the distribution
        virtual void setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance) {}

    protected:
        typeInt dim_;
        Vector mean_;
        Matrix cov_;
        Matrix covChol_;
        std::vector<PolynomialFamily> polyFamily_;
    };

    // Alias
    typedef std::shared_ptr<Distribution> DistributionPtr;
    typedef std::shared_ptr<const Distribution> DistributionConstPtr;

    DistributionPtr Dist(typeInt dim);
    DistributionPtr Dist(VectorConstRef mean, MatrixConstRef covariance, const std::vector<PolynomialFamily>& poly);
    DistributionPtr Dist(VectorConstRef mean, MatrixConstRef covariance, MatrixConstRef covChol,  const std::vector<PolynomialFamily>& poly);

    // Return the sum of dimensions of all elements of the distribution vector
    typeInt numberOfDimensions(const std::vector<DistributionPtr>& distributions);

}

#endif // DISTRIBUTION_HPP
