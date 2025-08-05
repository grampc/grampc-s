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


#include "distribution/univariate_gamma_distribution.hpp"

namespace grampc
{
    UnivariateGammaDistribution::UnivariateGammaDistribution(typeRNum alpha, typeRNum beta)
        : Distribution(Vector::Constant(1, 1, alpha / beta), 
                       Matrix::Constant(1, 1, alpha / (beta*beta)),
                       {PolynomialFamily::NONE}),
          dist_(alpha, beta),
          sample_(dim_)
    {
    }

    const Vector& UnivariateGammaDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }
    
    DistributionPtr Gamma(typeRNum alpha, typeRNum beta)
    {
        return DistributionPtr(new UnivariateGammaDistribution(alpha, beta));
    }
}