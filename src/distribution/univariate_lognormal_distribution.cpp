#include "distribution/univariate_lognormal_distribution.hpp"

namespace grampc
{
    UnivariateLognormalDistribution::UnivariateLognormalDistribution(typeRNum mu, typeRNum sigma)
        : Distribution(Vector::Constant(1, 1, std::exp(mu + 0.5 * sigma * sigma)), 
                       Matrix::Constant(1, 1, (std::exp(sigma * sigma) - 1.0) * std::exp(2.0 * mu + sigma * sigma)),
                       {PolynomialFamily::NONE}),
          dist_(mu, sigma),
          sample_(dim_)
    {
    }

    const Vector& UnivariateLognormalDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr LogNormal(typeRNum mu, typeRNum sigma)
    {
        return DistributionPtr(new UnivariateLognormalDistribution(mu, sigma));
    }
}