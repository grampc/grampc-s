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


#include "distribution/univariate_uniform_distribution.hpp"

namespace grampc
{
    UnivariateUniformDistribution::UnivariateUniformDistribution(typeRNum lowerBound, typeRNum upperBound)
        : Distribution(Vector::Constant(1, 1, (lowerBound + upperBound) * 0.5), 
                       Matrix::Constant(1, 1, (upperBound - lowerBound) * (upperBound - lowerBound) / 12.0),
                       {PolynomialFamily::LEGENDRE}),
          dist_(lowerBound, upperBound),
          sample_(dim_)
    {
    }

    const Vector& UnivariateUniformDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    void UnivariateUniformDistribution::setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance)
    {
        mean_ = mean;
        cov_ = covariance;
        covChol_ = cov_.llt().matrixL();

        // compute lower and upper bound
        typeRNum halfWidth = std::sqrt(12.0 * covariance(0,0)) * 0.5;
        dist_ = std::uniform_real_distribution<typeRNum>(mean(0) - halfWidth, mean(0) + halfWidth);
    }
   
    DistributionPtr Uniform(typeRNum lowerBound, typeRNum upperBound)
    {
        return DistributionPtr(new UnivariateUniformDistribution(lowerBound, upperBound));
    } 
}