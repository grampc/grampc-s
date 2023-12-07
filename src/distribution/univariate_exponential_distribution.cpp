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