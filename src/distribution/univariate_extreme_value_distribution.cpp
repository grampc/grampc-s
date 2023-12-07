#include "distribution/univariate_extreme_value_distribution.hpp"

namespace grampc
{
    UnivariateExtremeValueDistribution::UnivariateExtremeValueDistribution(typeRNum a, typeRNum b)
        : Distribution(Vector::Constant(1, 1, a + b * EULER_MASCHERONI_CONSTANT), 
                       Matrix::Constant(1, 1, PI * PI / 6.0 * b * b),
                       {PolynomialFamily::NONE}),
          dist_(a, b),
          sample_(dim_)
    {
    }

    const Vector& UnivariateExtremeValueDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr ExtremeValue(typeRNum a, typeRNum b)
    {
        return DistributionPtr(new UnivariateExtremeValueDistribution(a, b));
    }
}