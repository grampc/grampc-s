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


#ifndef UNIVARIATE_UNIFORM_DISTRIBUTION_HPP
#define UNIVARIATE_UNIFORM_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate uniform distribution p(x) = 1 / (upperBound - lowerBound) for lowerBound < x < upperBound
    class UnivariateUniformDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateUniformDistribution(typeRNum lowerBound, typeRNum upperBound);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;

        // Set new mean vector and covariance matrix of the distribution
        virtual void setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance) override;
    
    private:
        mutable std::uniform_real_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution
    DistributionPtr Uniform(typeRNum lowerBound, typeRNum upperBound);
}

#endif // UNIVARIATE_UNIFORM_DISTRIBUTION_HPP