#include "distribution/Gaussian_distribution.hpp"

namespace grampc
{
    GaussianDistribution::GaussianDistribution(VectorConstRef mean, MatrixConstRef covariance)
        : Distribution(mean, covariance, std::vector<PolynomialFamily>(mean.rows(), PolynomialFamily::HERMITE)),
          sample_(dim_),
          sampleNormalized_(dim_)
    {
    }

    GaussianDistribution::GaussianDistribution(typeRNum mean, typeRNum variance)
        : GaussianDistribution(Vector::Constant(1, 1, mean), Matrix::Constant(1, 1, variance))
    {
    }

    const Vector& GaussianDistribution::sample(RandomNumberGenerator &rng) const
    {
        // Get a sample from the multivariate standard normal distribution
        for (typeInt j = 0; j < dim_; ++j)
        {
            sampleNormalized_(j) = dist_(rng);
        }

        // Apply mean and covariance
        sample_ = mean_ + covChol_ * sampleNormalized_;
        return sample_;
    }

    void GaussianDistribution::setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance)
    {
        mean_ = mean;
        cov_ = covariance;
        covChol_ = cov_.llt().matrixL();
    }

    DistributionPtr Gaussian(VectorConstRef mean, MatrixConstRef covariance)
    {
        return DistributionPtr(new GaussianDistribution(mean, covariance));
    }

    DistributionPtr Gaussian(typeRNum mean, typeRNum variance)
    {
        return DistributionPtr(new GaussianDistribution(mean, variance));
    }
}