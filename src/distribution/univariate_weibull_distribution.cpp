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


#include "distribution/univariate_weibull_distribution.hpp"

namespace grampc
{
    UnivariateWeibullDistribution::UnivariateWeibullDistribution(typeRNum a, typeRNum b)
        : Distribution(Vector::Constant(1, 1, b * std::tgamma(1.0 + 1.0 / a)), 
                       Matrix::Constant(1, 1, b * b * (std::tgamma(1.0 + 2.0 / a) - std::tgamma(1.0 + 1.0 / a) * std::tgamma(1.0 + 1.0 / a))),
                       {PolynomialFamily::NONE}),
          dist_(a, b),
          sample_(dim_)
    {
    }

    const Vector& UnivariateWeibullDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr Weibull(typeRNum a, typeRNum b)
    {
        return DistributionPtr(new UnivariateWeibullDistribution(a, b));
    } 
}