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


#include "distribution/univariate_lognormal_distribution.hpp"

namespace grampc
{
    UnivariateLognormalDistribution::UnivariateLognormalDistribution(typeRNum mu, typeRNum sigma)
        : Distribution(Vector::Constant(1, 1, std::exp(mu + 0.5 * sigma * sigma)), 
                       Matrix::Constant(1, 1, (std::exp(sigma * sigma) - 1.0) * std::exp(2.0 * mu + sigma * sigma)),
                       {PolynomialFamily::NONE}),
          dist_(mu, sigma),
          sample_(dim_)
    {
    }

    const Vector& UnivariateLognormalDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr LogNormal(typeRNum mu, typeRNum sigma)
    {
        return DistributionPtr(new UnivariateLognormalDistribution(mu, sigma));
    }
}