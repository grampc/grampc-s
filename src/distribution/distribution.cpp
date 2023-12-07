#include "distribution/distribution.hpp"

namespace grampc
{
    Distribution::Distribution(typeInt dim)
        : dim_(dim),
          mean_(Vector::Zero(dim_)),
          cov_(Matrix::Zero(dim_, dim_)),
          covChol_(Matrix::Zero(dim_, dim_)),
          polyFamily_(dim, PolynomialFamily::NONE)
    {
    }

     Distribution::Distribution(VectorConstRef mean, MatrixConstRef covariance, const std::vector<PolynomialFamily>& poly)
        : dim_(mean.rows()),
          mean_(mean),
          cov_(covariance),
          covChol_(dim_, dim_),
          polyFamily_(poly)
    {
        // Compute the Cholesky decomposition of the covariance matrix
        covChol_ = cov_.llt().matrixL();
    }

    Distribution::Distribution(VectorConstRef mean, MatrixConstRef covariance, MatrixConstRef covChol,  const std::vector<PolynomialFamily>& poly)
        : dim_(mean.rows()),
          mean_(mean),
          cov_(covariance),
          covChol_(covChol),
          polyFamily_(poly)
    {
    }

    const Vector& Distribution::mean() const
    {
        return mean_;
    }

    const Matrix& Distribution::covariance() const
    {
        return cov_;
    }

    const Matrix& Distribution::covCholesky() const
    {
        return covChol_;
    }

    const Vector& Distribution::sample(RandomNumberGenerator& rng) const
    {
        std::cerr << "Sampling of an abstract distribution is not possible!" << std::endl;
        return mean_;
    }

    typeInt Distribution::dimension() const
    {
        return dim_;
    }

    const std::vector<PolynomialFamily> &Distribution::polynomialFamily() const
    {
        return polyFamily_;
    }

    DistributionPtr Dist(typeInt dim)
    {
        return DistributionPtr(new Distribution(dim));
    }

    typeInt numberOfDimensions(const std::vector<DistributionPtr>& distributions)
    {
        typeInt dim = 0;
        for(DistributionConstPtr dist : distributions)
        {
            dim += dist->dimension();
        }
        
        return dim;
    }
}