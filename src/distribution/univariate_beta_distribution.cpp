#include "distribution/univariate_beta_distribution.hpp"

namespace grampc
{
    UnivariateBetaDistribution::UnivariateBetaDistribution(typeRNum p, typeRNum q)
        : Distribution(Vector::Constant(1, 1, p / (p + q)), 
                       Matrix::Constant(1, 1, p * q / ((p + q + 1) * (p + q) * (p + q))),
                       {PolynomialFamily::NONE}),
          distGamma1_(p, 1),
          distGamma2_(q, 1),
          sample_(dim_)
    {
    }

    const Vector& UnivariateBetaDistribution::sample(RandomNumberGenerator& rng) const
    {
        temp_ = distGamma1_(rng);
        sample_[0] = temp_ / (temp_ + distGamma2_(rng));
        return sample_;
    }
    
    DistributionPtr Beta(typeRNum p, typeRNum q)
    {
        return DistributionPtr(new UnivariateBetaDistribution(p, q));
    }
}