#include "distribution/univariate_f_distribution.hpp"

namespace grampc
{
    UnivariateFDistribution::UnivariateFDistribution(typeRNum m, typeRNum n)
        : Distribution(Vector::Constant(1, 1, n / (n - 2.0)), 
                       Matrix::Constant(1, 1, (2.0 * n * n * (m + n - 2.0)) / (m * (n-2.0)*(n-2.0) * (n-4.0))),
                       {PolynomialFamily::NONE}),
          dist_(m, n),
          sample_(dim_)
    {
    }

    const Vector& UnivariateFDistribution::sample(RandomNumberGenerator& rng) const
    {
        sample_[0] = dist_(rng);
        return sample_;
    }

    DistributionPtr Fisher(typeRNum m, typeRNum n)
    {
        return DistributionPtr(new UnivariateFDistribution(m, n));
    }  
}