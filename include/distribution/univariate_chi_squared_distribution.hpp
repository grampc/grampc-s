#ifndef UNIVARIATE_CHI_SQUARED_DISTRIBUTION_HPP
#define UNIVARIATE_CHI_SQUARED_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate chi-squared distribution p(x) = (x^(n/2 - 1) * exp(-x/2)) / (gamma(n/2) * 2^(n/2)) with n > 0, x > 0
    class UnivariateChiSquaredDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateChiSquaredDistribution(typeRNum n);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::chi_squared_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    DistributionPtr ChiSquared(typeRNum n);  
}

#endif // UNIVARIATE_CHI_SQUARED_DISTRIBUTION_HPP