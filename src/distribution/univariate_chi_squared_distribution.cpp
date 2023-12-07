#include "distribution/univariate_chi_squared_distribution.hpp"

namespace grampc
{
    UnivariateChiSquaredDistribution::UnivariateChiSquaredDistribution(typeRNum n)
        : Distribution(Vector::Constant(1, 1, n), 
                       Matrix::Constant(1, 1, 2.0 * n),
                       {PolynomialFamily::NONE}),
          dist_(n),
          sample_(dim_)
    {
    }

    const Vector& UnivariateChiSquaredDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr ChiSquared(typeRNum n)
    {
        return DistributionPtr(new UnivariateChiSquaredDistribution(n));
    } 
}