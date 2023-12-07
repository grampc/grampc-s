#ifndef UNIVARIATE_GAMMA_DISTRIBUTION_HPP
#define UNIVARIATE_GAMMA_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate gamma distribution p(x) = exp(-x/beta) / (beta^alpha * gamma(alpha)) * x^(alpha - 1) with x > 0, alpha > 0, beta > 0
    class UnivariateGammaDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateGammaDistribution(typeRNum alpha, typeRNum beta);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::gamma_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    DistributionPtr Gamma(typeRNum alpha, typeRNum beta);
}

#endif // UNIVARIATE_GAMMA_DISTRIBUTION_HPP