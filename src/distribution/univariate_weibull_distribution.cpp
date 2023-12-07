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