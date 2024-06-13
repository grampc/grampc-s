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


#ifndef GAUSSIAN_DISTRIBUTION_HPP
#define GAUSSIAN_DISTRIBUTION_HPP

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

    private:
        mutable std::normal_distribution<typeRNum> dist_{0.0, 1.0};
        mutable Vector sample_;
        mutable Vector sampleNormalized_;
    };

    // Constructor for a multivariate Gaussian distribution
    DistributionPtr Gaussian(VectorConstRef mean, MatrixConstRef covariance);
    // Constructor for a univariate Gaussian distribution
    DistributionPtr Gaussian(typeRNum mean, typeRNum variance);
}

#endif // GAUSSIAN_DISTRIBUTION_HPP