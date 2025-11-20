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


#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    MultivariateDistribution::MultivariateDistribution(const std::vector<DistributionPtr>& distributions)
        : Distribution(numberOfDimensions(distributions)),
          numDistributions_((typeInt) distributions.size()),
          distributions_(distributions),
          dimensionVec_(numDistributions_),
          sample_(dim_)
    {
        // Counters for dimensions and distributions
        typeInt dimCounter = 0;
        typeInt distCounter = 0;

        // Iterator for polynomial families
        std::vector<PolynomialFamily>::iterator polyIter = polyFamily_.begin();

        // Go through all distributions
        for(DistributionPtr& dist : distributions_)
        {
            // Dimension of the current distribution
            dimensionVec_[distCounter] = dist->dimension();

            //Add mean
            mean_.segment(dimCounter, dimensionVec_[distCounter]) = dist->mean();

            // Add covariance matrix and its Cholesky decomposition
            cov_.block(dimCounter, dimCounter, dimensionVec_[distCounter], dimensionVec_[distCounter]) = dist->covariance();
            covChol_.block(dimCounter, dimCounter, dimensionVec_[distCounter], dimensionVec_[distCounter]) = dist->covCholesky();

            // Add polynomial family
            const std::vector<PolynomialFamily>& tempPoly = dist->polynomialFamily();
            polyIter = std::copy(tempPoly.begin(), tempPoly.end(), polyIter);

            // Increase counters
            dimCounter += dimensionVec_[distCounter];
            distCounter++;
        }
    }

    const Vector& MultivariateDistribution::sample(RandomNumberGenerator& rng) const
    {
        typeInt iterator = 0;

         // Go through all distributions
        for(const DistributionConstPtr& dist : distributions_)
        {
            if(dist->dimension() > 0)
            {
                sample_.segment(iterator, dist->dimension()) = dist->sample(rng);
                iterator += dist->dimension();
            }
        }
        return sample_;
    }

    void MultivariateDistribution::setMeanAndCovariance(VectorConstRef mean, MatrixConstRef covariance)
    {
        // Counters for dimensions and distributions
        typeInt dimCounter = 0;
        typeInt distCounter = 0;

        // Go through all distributions
        for(DistributionPtr& dist : distributions_)
        {
            // Dimension of the current distribution
            dimensionVec_[distCounter] = dist->dimension();

            // set mean and covariance
            dist->setMeanAndCovariance(mean.segment(dimCounter, dimensionVec_[distCounter]), covariance.block(dimCounter, dimCounter, dimensionVec_[distCounter], dimensionVec_[distCounter]));

            //Add mean
            mean_.segment(dimCounter, dimensionVec_[distCounter]) = dist->mean();

            // Add covariance matrix and its Cholesky decomposition
            cov_.block(dimCounter, dimCounter, dimensionVec_[distCounter], dimensionVec_[distCounter]) = dist->covariance();
            covChol_.block(dimCounter, dimCounter, dimensionVec_[distCounter], dimensionVec_[distCounter]) = dist->covCholesky();

            // Increase counters
            dimCounter += dimensionVec_[distCounter];
            distCounter++;
        }
    }


    void MultivariateDistribution::replaceDistribution(typeInt index, DistributionPtr newDist)
    {
        typeInt distDim = dimensionVec_[index];
        typeInt dimCounter = dimensionVec_.head(index).sum();

        // Check dimensions
        if(newDist->dimension() != distDim)
        {
            std::cerr << "Dimension of the new distribution does not match dimension of the current distribution!" << std::endl;
        }

        // Replace distribution
        distributions_[index] = newDist;

        // Replace vectors and matrices
        mean_.segment(dimCounter, distDim) = newDist->mean();
        cov_.block(dimCounter, dimCounter, distDim, distDim) = newDist->covariance();
        covChol_.block(dimCounter, dimCounter, distDim, distDim) = newDist->covCholesky();

        // Replace polynomial family
        const std::vector<PolynomialFamily>& tempPoly = newDist->polynomialFamily();
        std::copy(tempPoly.begin(), tempPoly.end(), polyFamily_.begin() + dimCounter);
    }

    DistributionConstPtr MultivariateDistribution::getDistribution(typeInt index) const
    {
        return distributions_[index];
    }

    MultiDistributionPtr MultiDist(const std::vector<DistributionPtr>& distributions)
    {
        return MultiDistributionPtr(new MultivariateDistribution(distributions));
    } 
}