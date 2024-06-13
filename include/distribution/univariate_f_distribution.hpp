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


#ifndef UNIVARIATE_F_DISTRIBUTION_HPP
#define UNIVARIATE_F_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate chi-squared distribution p(x) = gamma((m+n)/2) / (gamma(m/2) * gamma(n/2)) * (m/n)^(m/2) * x^(m/2 - 1) * (1 + m/n * x)^(-(m+n)/2) with x > 0, m > 0, n > 4
    class UnivariateFDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution
        UnivariateFDistribution(typeRNum m, typeRNum n);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::fisher_f_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution
    DistributionPtr Fisher(typeRNum m, typeRNum n); 
}

#endif // UNIVARIATE_CHI_SQUARED_DISTRIBUTION_HPP