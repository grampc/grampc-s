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


#include "distribution/univariate_student-t_distribution.hpp"

namespace grampc
{
    UnivariateStudentTDistribution::UnivariateStudentTDistribution(typeRNum nu, typeRNum mu, typeRNum sigma)
        : Distribution(Vector::Constant(1, 1, mu), 
                       Matrix::Constant(1, 1, sigma * sigma * nu / (nu - 2.0)),
                       {PolynomialFamily::NONE}),
          dist_(nu),
          sample_(dim_),
          mu_(mu),
          sigma_(sigma)
    {
    }

    const Vector& UnivariateStudentTDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = mu_ + sigma_ * dist_(rng);
        return sample_;
    }

    DistributionPtr StudentT(typeRNum nu)
    {
        return DistributionPtr(new UnivariateStudentTDistribution(nu));
    } 

    DistributionPtr StudentT(typeRNum nu, typeRNum mu, typeRNum sigma)
    {
        return DistributionPtr(new UnivariateStudentTDistribution(nu, mu, sigma));
    }
}