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


#ifndef UNIVARIATE_PIECEWISE_CONSTANT_DISTRIBUTION_HPP
#define UNIVARIATE_PIECEWISE_CONSTANT_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate piecewise constant distribution. The intervalLimits vector starts with the beginning of the first interval and ends with the end of the last interval
    class UnivariatePiecewiseConstantDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariatePiecewiseConstantDistribution(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        typeRNum mean(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity);
        typeRNum variance(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity);

        std::vector<typeRNum> intervalLimits_;
        std::vector<typeRNum> intervalProbabilityDensity_;
        std::vector<typeRNum> weights_;
        mutable std::piecewise_constant_distribution<typeRNum> dist_;
        mutable Vector sample_;
        typeInt numIntervals_;
    };
    
    // Constructor for the univariate distribution 
    DistributionPtr PiecewiseConstant(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity);
}

#endif // UNIVARIATE_PIECEWISE_CONSTANT_DISTRIBUTION_HPP