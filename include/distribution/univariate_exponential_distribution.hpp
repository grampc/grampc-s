#ifndef UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP
#define UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate exponential distribution p(x) = lambda * exp(-lambda * x) with lambda > 0, x >= 0
    class UnivariateExponentialDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateExponentialDistribution(typeRNum lambda);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::exponential_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    DistributionPtr Exponential(typeRNum lambda);
}

#endif // UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP