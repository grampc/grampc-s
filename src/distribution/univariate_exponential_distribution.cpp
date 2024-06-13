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


#include "distribution/univariate_exponential_distribution.hpp"

namespace grampc
{
    UnivariateExponentialDistribution::UnivariateExponentialDistribution(typeRNum lambda)
        : Distribution(Vector::Constant(1, 1, 1.0 / lambda), 
                       Matrix::Constant(1, 1, 1.0 / (lambda * lambda)),
                       {PolynomialFamily::NONE}),
          dist_(lambda),
          sample_(dim_)
    {
    }

    const Vector& UnivariateExponentialDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr Exponential(typeRNum lambda)
    {
        return DistributionPtr(new UnivariateExponentialDistribution(lambda));
    }  
}