#include "point_transformation/monte_carlo.hpp"


namespace grampc
{
    MonteCarloTransformation::MonteCarloTransformation(typeInt dimX, typeInt dimY, typeInt numberOfPoints, const RandomNumberGenerator& rng)
    : dimX_(dimX),
      dimY_(dimY),
      rng_(rng),
      numberOfPoints_(numberOfPoints),
      points_(dimX, numberOfPoints),
      meanX_(dimX),
      meanY_(dimY),
      covariance_(dimX, dimY),
      dmean1D_dpoints_(numberOfPoints),
      dvar_dpoints_(numberOfPoints_),
      weight_(1.0 / ((typeRNum) numberOfPoints)),
      diffToMeanX_(dimX, numberOfPoints),
      diffToMeanY_(dimY, numberOfPoints),
      diffToMean1D_(numberOfPoints)
    { 
        dmean1D_dpoints_.setConstant(weight_);
    }

    const Matrix& MonteCarloTransformation::points(DistributionConstPtr dist)
    {
        for(typeInt pointIndex = 0; pointIndex < numberOfPoints_; ++pointIndex)
        {
            // Samples random points from distribution
            points_.col(pointIndex) = dist->sample(rng_);
        }
        return points_;
    }

    const Matrix& MonteCarloTransformation::points()
    {
        return points_;
    }

    const Vector& MonteCarloTransformation::mean(MatrixConstRef points)
    {
        meanY_ = points.rowwise().mean();
        return meanY_;
    }

    typeRNum MonteCarloTransformation::mean1D(RowVectorConstRef points)
    {
        return points.mean();
    }

    const Matrix& MonteCarloTransformation::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        // Compute mean of points
        meanX_.noalias() = pointsX.rowwise().mean();
        meanY_.noalias() = pointsY.rowwise().mean();

        // difference between points and mean
        for(typeInt i = 0; i < numberOfPoints_; ++i)
        {
            diffToMeanX_.col(i) = (pointsX.col(i) - meanX_) * weight_;
            diffToMeanY_.col(i) = pointsY.col(i) - meanY_;
        }

        // Compute and return the covariance matrix
        covariance_.noalias() = diffToMeanX_ * diffToMeanY_.transpose();
        return covariance_;
    }

    typeRNum MonteCarloTransformation::variance(RowVectorConstRef points)
    {
        tempScalar = points.mean();
        for(typeInt i = 0; i < numberOfPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        return weight_ * diffToMean1D_ * diffToMean1D_.transpose();
    }

    const Vector& MonteCarloTransformation::dmean1D_dpoints()
    {
        return dmean1D_dpoints_;
    }

    const Vector& MonteCarloTransformation::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points.mean();

        // difference between points and mean
        for(typeInt i = 0; i < numberOfPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }

        tempScalar = - 2.0 * weight_ * diffToMean1D_.mean(); 
        
        // Derivative of the variance
        for(typeInt i = 0; i < numberOfPoints_; ++i)
        {
            dvar_dpoints_(i) = 2.0 * weight_ * diffToMean1D_(i) + tempScalar;
        }
        return dvar_dpoints_;
    }

    typeInt MonteCarloTransformation::numberOfPoints() const
    {
        return numberOfPoints_;
    }

    PointTransformationPtr MonteCarlo(typeInt dimX, typeInt dimY, typeInt numberOfPoints, const RandomNumberGenerator& rng)
    {
        return PointTransformationPtr(new MonteCarloTransformation(dimX, dimY, numberOfPoints, rng));
    }
}