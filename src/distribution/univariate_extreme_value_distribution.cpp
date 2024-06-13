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


#include "distribution/univariate_extreme_value_distribution.hpp"

namespace grampc
{
    UnivariateExtremeValueDistribution::UnivariateExtremeValueDistribution(typeRNum a, typeRNum b)
        : Distribution(Vector::Constant(1, 1, a + b * EULER_MASCHERONI_CONSTANT), 
                       Matrix::Constant(1, 1, PI * PI / 6.0 * b * b),
                       {PolynomialFamily::NONE}),
          dist_(a, b),
          sample_(dim_)
    {
    }

    const Vector& UnivariateExtremeValueDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr ExtremeValue(typeRNum a, typeRNum b)
    {
        return DistributionPtr(new UnivariateExtremeValueDistribution(a, b));
    }
}