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


#ifndef UNIVARIATE_EXTREME_VALUE_DISTRIBUTION_HPP
#define UNIVARIATE_EXTREME_VALUE_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"
#include "util/grampc_s_constants.hpp"

namespace grampc
{
    // Univariate extreme value distribution p(x) = 1 / b * exp((a − x) / b − exp((a − x) / b))
    class UnivariateExtremeValueDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateExtremeValueDistribution(typeRNum a, typeRNum b);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::extreme_value_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution 
    DistributionPtr ExtremeValue(typeRNum a, typeRNum b);
}

#endif // UNIVARIATE_EXTREME_VALUE_DISTRIBUTION_HPP