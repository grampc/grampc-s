#ifndef MULTIVARIATE_UNCORRELATED_DISTRIBUTION_HPP
#define MULTIVARIATE_UNCORRELATED_DISTRIBUTION_HPP

#include "distribution.hpp"
#include "util/grampc_s_util.hpp"

namespace grampc
{
    // Multivariate distribution consisting of multiple (univariate or multivariate) uncorrelated distributions p(x) = p_1(x_1) * p_2(x_2) * p_3(x_3) ...
    class MultivariateDistribution : public Distribution
    {
    public:
        // Constructor using uncorrelated distributions
        MultivariateDistribution(const std::vector<DistributionPtr>& distributions);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;

        // Set new mean vector and covariance matrix of the distribution
        virtual void setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance) override;

        // Replace the distribution with specified index by a new distribution with same dimension
        void replaceDistribution(typeInt index, DistributionPtr newDist);

        DistributionConstPtr getDistribution(typeInt index) const;
        
    private:
        typeInt numDistributions_;
        std::vector<DistributionPtr> distributions_;
        Vector dimensionVec_;
        mutable Vector sample_;
    };

    // Alias
    typedef std::shared_ptr<MultivariateDistribution> MultiDistributionPtr;
    typedef std::shared_ptr<const MultivariateDistribution> MultiDistributionConstPtr;

    MultiDistributionPtr MultiDist(const std::vector<DistributionPtr>& distributions);
}

#endif // MULTIVARIATE_UNCORRELATED_DISTRIBUTION_HPP