#ifndef GAUSSIAN_DISTRIBUTION_HPP
#define GAUSSIAN_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"
#include "util/grampc_s_util.hpp"

namespace grampc
{
    // Multivariate Gaussian distribution p(x) = (2*pi)^(-d/2) * det(Sigma)^(-1/2) * exp(-1/2 * (x - mu)^T * Sigma^(-1) * (x - mu))
    // with mean vector mu, covariance matrix Sigma, and dimension d of the random variable
    class GaussianDistribution : public Distribution
    {
    public:
        // Constructor for a multivariate Gaussian distribution
        GaussianDistribution(VectorConstRef mean, MatrixConstRef covariance);

        // Constructor for a univariate Gaussian distribution
        GaussianDistribution(typeRNum mean, typeRNum variance);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator &rng) const override;

        // Set new mean vector and covariance matrix of the distribution
        virtual void setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance) override;

    private:
        mutable std::normal_distribution<typeRNum> dist_{0.0, 1.0};
        mutable Vector sample_;
        mutable Vector sampleNormalized_;
    };

    DistributionPtr Gaussian(VectorConstRef mean, MatrixConstRef covariance);
    DistributionPtr Gaussian(typeRNum mean, typeRNum variance);
}

#endif // GAUSSIAN_DISTRIBUTION_HPP