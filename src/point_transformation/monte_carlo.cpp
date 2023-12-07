#include "point_transformation/monte_carlo.hpp"


namespace grampc
{
    MonteCarloTransformation::MonteCarloTransformation(typeInt dim, typeInt numberOfPoints, const RandomNumberGenerator& rng)
    : rng_(rng),
      numberOfPoints_(numberOfPoints),
      dim_(dim),
      points_(dim, numberOfPoints),
      mean_(dim),
      covariance_(dim, dim),
      dmean_dpoints_vec_(dim * numberOfPoints),
      dmean1D_dpoints_(numberOfPoints),
      dcov_dpoints_vec_(dim * numberOfPoints),
      dvar_dpoints_(numberOfPoints_),
      weight_(1.0 / ((typeRNum) numberOfPoints)),
      diffToMean_(dim, numberOfPoints),
      diffToMean1D_(numberOfPoints),
      dcov_(dim, numberOfPoints),
      diffToMean_jk_dx_(Matrix::Zero(dim, numberOfPoints)),
      diffToMean_ik_dx_(Matrix::Zero(dim, numberOfPoints))
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

    const Vector& MonteCarloTransformation::mean(MatrixConstRef points)
    {
        mean_ = points.rowwise().mean();
        return mean_;
    }

    typeRNum MonteCarloTransformation::mean1D(RowVectorConstRef points)
    {
        return points.mean();
    }

    const Matrix& MonteCarloTransformation::covariance(MatrixConstRef points)
    {
        // Compute the mean
        mean_ = points.rowwise().mean();

        // Difference between points and mean
        diffToMean_ = points.colwise() - mean_;

        // Compute and return the covariance matrix
        covariance_.noalias() = weight_ * diffToMean_ * diffToMean_.transpose();
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

    const Vector& MonteCarloTransformation::dmean_dpoints_vec(VectorConstRef vec)
    {
        // out = (w_1*v_1  w_1*v_2  w_1*v_3  ...  w_2*v_1  w_2*v_2   w_2*v_3  ...)^T

        // Set and return output
        dmean_dpoints_vec_ = vec.replicate(numberOfPoints_, 1) * weight_;
        return dmean_dpoints_vec_;
    }

    const Vector& MonteCarloTransformation::dmean1D_dpoints()
    {
        return dmean1D_dpoints_;
    }

    const Vector& MonteCarloTransformation::dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec)
    {
        /*
        * Compute out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        * with dcov_ij/dx = sum_k(weightsVar_k * (x_ik - m_i) * (dx_jk/dx - dm_j/d_x) +
        *                         weightsVar_k * (x_jk - m_j) * (dx_ik/dx - dm_i/d_x))
        */

        // Mapping of the vec input to a matrix corresponding to the covariance matrix
        Eigen::Map<const Matrix> vecMatrix(vec.data(), dim_, dim_);

        // Compute the mean and the difference from each sigma point to the mean
        mean_ = points.rowwise().mean();
        diffToMean_ = points.colwise() - mean_;

        // out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        // Initialize output
        dcov_dpoints_vec_.setZero();             
        
        // Compute dcov_ij/dx
        for(typeInt j = 0; j < dim_; ++j)
        {
            // diffToMean_jk_dx = "1 in element jk" + "weightsMean in row j"
            diffToMean_jk_dx_.row(j).setConstant(weight_);

            for(typeInt i = 0; i < dim_; ++i)
            {
                // diffToMean_ik_dx = "1 in element ik" + "weightsMean in row i"
                diffToMean_ik_dx_.row(i).setConstant(weight_);
                
                // Set dcov to zero
                dcov_.setZero();

                // Go through all sigma points
                for(typeInt k = 0; k < numberOfPoints_; ++k)
                {
                    // Add 1 at elements jk and ik
                    diffToMean_jk_dx_(j, k) += 1.0;
                    diffToMean_ik_dx_(i, k) += 1.0;

                    // Add dcov_ij_dx for sigmapoint k
                    dcov_ += weight_ * (diffToMean_(i, k) * diffToMean_jk_dx_ + diffToMean_(j, k) * diffToMean_ik_dx_);

                    // Remove 1 at elements jk and ik
                    diffToMean_jk_dx_(j, k) -= 1.0;
                    diffToMean_ik_dx_(i, k) -= 1.0;
                }
                // Remove row j
                diffToMean_jk_dx_.row(i).setZero();

                // Add vecMatrix_ij * dcov_ij/dx to the output
                dcov_dpoints_vec_ += vecMatrix(i, j) * dcov_.reshaped(dim_ * numberOfPoints_, 1);
            }
            // Remove row i
            diffToMean_ik_dx_.row(j).setZero();
        }
        return dcov_dpoints_vec_;
    }

    const Vector& MonteCarloTransformation::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points.mean();

        // Difference of each point to the mean
        for(typeInt i = 0; i < numberOfPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        
        // derivative of the variance
        tempScalar = diffToMean1D_.sum() * weight_;   
        dvar_dpoints_ = 2.0 * weight_ * (diffToMean1D_.array() - tempScalar).matrix();
        return dvar_dpoints_;
    }

    typeInt MonteCarloTransformation::numberOfPoints() const
    {
        return numberOfPoints_;
    }

    PointTransformationPtr MonteCarlo(typeInt dim, typeInt numberOfPoints, const RandomNumberGenerator& rng)
    {
        return PointTransformationPtr(new MonteCarloTransformation(dim, numberOfPoints, rng));
    }
}