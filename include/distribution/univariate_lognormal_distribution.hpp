#ifndef UNIVARIATE_LOGNORMAL_DISTRIBUTION_HPP
#define UNIVARIATE_LOGNORMAL_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate log-normal distribution p(x) = 1/(sigma * x * sqrt(2*pi)) * exp(−(ln(x)-mu)^2 / (2*sigma^2)) for x > 0
    class UnivariateLognormalDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateLognormalDistribution(typeRNum mu, typeRNum sigma);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::lognormal_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    DistributionPtr LogNormal(typeRNum mu, typeRNum sigma);
}

#endif // UNIVARIATE_LOGNORMAL_DISTRIBUTION_HPP