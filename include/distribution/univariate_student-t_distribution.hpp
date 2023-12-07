#ifndef UNIVARIATE_STUDENT_T_DISTRIBUTION_HPP
#define UNIVARIATE_STUDENT_T_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"


namespace grampc
{
    // Univariate student-t distribution p(x) = 1/(sigma * sqrt(nu*pi)) * gamma((nu+1)/2) / gamma(nu/2) * (1 + ((x-mu)/sigma)^2 / nu)^(−(nu+1)/2) with nu > 2, location parameter mu, and scale parameter sigma
    class UnivariateStudentTDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution
        UnivariateStudentTDistribution(typeRNum nu, typeRNum mu = 0.0, typeRNum sigma = 1.0);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::student_t_distribution<typeRNum> dist_;
        mutable Vector sample_;
        typeRNum mu_;
        typeRNum sigma_;
    };
    
    DistributionPtr StudentT(typeRNum nu);
    DistributionPtr StudentT(typeRNum nu, typeRNum mu, typeRNum sigma);    
}

#endif // UNIVARIATE_STUDENT_T_DISTRIBUTION_HPP