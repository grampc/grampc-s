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


#include "distribution/univariate_piecewise_constant_distribution.hpp"

namespace grampc
{
    UnivariatePiecewiseConstantDistribution::UnivariatePiecewiseConstantDistribution(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity)
        : Distribution(Vector::Constant(1, 1, mean(intervalLimits, intervalProbabilityDensity)), 
                       Matrix::Constant(1, 1, variance(intervalLimits, intervalProbabilityDensity)),
                       {PolynomialFamily::NONE}),
          intervalLimits_(intervalLimits),
          intervalProbabilityDensity_(intervalProbabilityDensity),
          weights_(intervalProbabilityDensity.size()),
          sample_(dim_),
          numIntervals_(intervalProbabilityDensity_.size())
    {
        // Compute weights
        for(typeInt i = 0; i < intervalProbabilityDensity_.size(); ++i)
        {
            weights_[i] = intervalProbabilityDensity_[i] * (intervalLimits_[i+1] - intervalLimits_[i]);
        }

        // Set distribution
        dist_ = std::piecewise_constant_distribution<typeRNum>(intervalLimits_.begin(), intervalLimits_.end(), weights_.begin());
    }

    const Vector& UnivariatePiecewiseConstantDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    typeRNum UnivariatePiecewiseConstantDistribution::mean(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity)
    {
        typeRNum out = 0.0;
        typeInt numIntervals = intervalProbabilityDensity.size();

        // Compute the mean of the distribution as E{X} = sum(w_i * E{U_i}), where U_i is the uniform distribution of one interval and w_i is its probability
        for(typeInt i = 0; i < numIntervals; ++i)
        {
            out += intervalProbabilityDensity[i] * (intervalLimits[i+1] - intervalLimits[i]) * (intervalLimits[i+1] + intervalLimits[i]) * 0.5;
        }
        return out;
    }

    typeRNum UnivariatePiecewiseConstantDistribution::variance(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity)
    {
        /************************************************************************************************************************************************
        * Compute the variance of the distribution as Var{X} = E{X^2} - E{X}^2 
        * with E{X^2} = sum(w_i * E{U_i^2})
        * For each interval with boundaries a and b: E{U_i^2} = Var{U_i} + E{U_i}^2 = 1/12 * (b-a)^2 + ((b+a)/2)^2 = 1/12 * (b-a)^2 + (b+a)^2 / 4
        ************************************************************************************************************************************************/

        typeInt numIntervals = intervalProbabilityDensity.size();

        // Vector of second moments E{U_i^2} of uniform distributions with boundaries a and b
        std::vector<typeRNum> secondMoments(numIntervals);

        // Boundaries
        typeRNum a;
        typeRNum b;

        // Second central moment of the piecewise constant distribution E{X^2}
        typeRNum secondMomentTotal = 0.0;

        // Compute the second moment E{U_i^2} = 1/12 * (b-a)^2 + (b+a)^2 / 4 for each uniform distribution
        for(typeInt i = 0; i < numIntervals; ++i)
        {
            a = intervalLimits[i];
            b = intervalLimits[i+1];
            secondMoments[i] = ((b-a)*(b-a) / 12.0) + (b+a)*(b+a) / 4.0;
        }

        // Compute E{X^2} = sum(w_i * E{U_i^2})
        for(typeInt i = 0; i < numIntervals; ++i)
        {
            secondMomentTotal += intervalProbabilityDensity[i] * (intervalLimits[i+1] - intervalLimits[i]) * secondMoments[i];
        }

        // Compute Variance of the piecewise constant distribution Var{X} = E{X^2} - E{X}^2
        typeRNum m = this->mean(intervalLimits, intervalProbabilityDensity);
        return secondMomentTotal - m * m;
    }

    DistributionPtr PiecewiseConstant(const std::vector<typeRNum>& intervalLimits, const std::vector<typeRNum>& intervalProbabilityDensity)
    {
        return DistributionPtr(new UnivariatePiecewiseConstantDistribution(intervalLimits, intervalProbabilityDensity));
    }
}