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


#ifndef UNIVARIATE_WEIBULL_DISTRIBUTION_HPP
#define UNIVARIATE_WEIBULL_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate Weibull distribution p(x) = a/b * (x/b)^(a−1) * exp(−(x/b)^a) for x > 0
    class UnivariateWeibullDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution
        UnivariateWeibullDistribution(typeRNum a, typeRNum b);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::weibull_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution
    DistributionPtr Weibull(typeRNum a, typeRNum b);
}

#endif // UNIVARIATE_WEIBULL_DISTRIBUTION_HPP