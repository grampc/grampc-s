#ifndef UNIVARIATE_BETA_DISTRIBUTION_HPP
#define UNIVARIATE_BETA_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate beta distribution p(x) = 1 / B(p, q) * x^(p-1) * (1 - x)^(q - 1) with 0 < x < 1, p > 0, q > 0
    class UnivariateBetaDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateBetaDistribution(typeRNum p, typeRNum q);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable typeRNum temp_;
        mutable std::gamma_distribution<typeRNum> distGamma1_;
        mutable std::gamma_distribution<typeRNum> distGamma2_;
        mutable Vector sample_;
    };
    
    DistributionPtr Beta(typeRNum p, typeRNum q);
}

#endif // UNIVARIATE_BETA_DISTRIBUTION_HPP