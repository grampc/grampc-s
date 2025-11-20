/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


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
        : dim_((typeInt) mean.rows()),
          mean_(mean),
          cov_(covariance),
          covChol_(dim_, dim_),
          polyFamily_(poly)
    {
        // Compute the Cholesky decomposition of the covariance matrix
        covChol_ = cov_.llt().matrixL();
    }

    Distribution::Distribution(VectorConstRef mean, MatrixConstRef covariance, MatrixConstRef covChol,  const std::vector<PolynomialFamily>& poly)
        : dim_((typeInt) mean.rows()),
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

    void Distribution::setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance)
    {
        mean_ = mean;
        cov_ = covariance;
        covChol_ = cov_.llt().matrixL();
    }

    DistributionPtr Dist(typeInt dim)
    {
        return DistributionPtr(new Distribution(dim));
    }

    DistributionPtr Dist(VectorConstRef mean, MatrixConstRef covariance, const std::vector<PolynomialFamily>& poly)
    {
        return DistributionPtr(new Distribution(mean, covariance, poly));
    }

    DistributionPtr Dist(VectorConstRef mean, MatrixConstRef covariance, MatrixConstRef covChol,  const std::vector<PolynomialFamily>& poly)
    {
        return DistributionPtr(new Distribution(mean, covariance, covChol, poly));
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