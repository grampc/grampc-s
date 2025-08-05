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


#include "distribution/univariate_chi_squared_distribution.hpp"

namespace grampc
{
    UnivariateChiSquaredDistribution::UnivariateChiSquaredDistribution(typeRNum n)
        : Distribution(Vector::Constant(1, 1, n), 
                       Matrix::Constant(1, 1, 2.0 * n),
                       {PolynomialFamily::NONE}),
          dist_(n),
          sample_(dim_)
    {
    }

    const Vector& UnivariateChiSquaredDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr ChiSquared(typeRNum n)
    {
        return DistributionPtr(new UnivariateChiSquaredDistribution(n));
    } 
}