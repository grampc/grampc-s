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


#include "distribution/univariate_f_distribution.hpp"

namespace grampc
{
    UnivariateFDistribution::UnivariateFDistribution(typeRNum m, typeRNum n)
        : Distribution(Vector::Constant(1, 1, n / (n - 2.0)), 
                       Matrix::Constant(1, 1, (2.0 * n * n * (m + n - 2.0)) / (m * (n-2.0)*(n-2.0) * (n-4.0))),
                       {PolynomialFamily::NONE}),
          dist_(m, n),
          sample_(dim_)
    {
    }

    const Vector& UnivariateFDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr Fisher(typeRNum m, typeRNum n)
    {
        return DistributionPtr(new UnivariateFDistribution(m, n));
    }  
}