#include "point_transformation/unscented_transformation.hpp"

namespace grampc
{
    UnscentedTransformation::UnscentedTransformation(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain)
    : alpha_(alpha),
      beta_(beta),
      kappa_(kappa),
      dim_(dim),
      considerUncertain_(considerUncertain),
      numUncertainVariables_(std::count(considerUncertain_.begin(), considerUncertain_.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dim, numPoints_),
      points_(dim, numPoints_),
      mean_(dim),
      covariance_(dim, dim),
      dmean_dpoints_vec_(dim*numPoints_),
      dcov_dpoints_vec_(dim*numPoints_),
      dvar_dpoints_(numPoints_),
      weightsMean_(numPoints_),
      weightsVar_(numPoints_),
      diffToMean_(dim, numPoints_),
      diffToMean1D_(numPoints_),
      dcov_(dim, numPoints_),
      diffToMean_jk_dx_(Matrix::Zero(dim, numPoints_)),
      diffToMean_ik_dx_(Matrix::Zero(dim, numPoints_))
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

        // Index for uncertain variables
        typeInt uncertIndex = 0;

        // Fill vector with indices of uncertain variables
        for(typeInt i = 0; i < dim; ++i)
        {
            if(considerUncertain_[i])
            {
                uncertainIndices[uncertIndex++] = i;
            }
        }

        // Set normalized Points
        Matrix A = Matrix::Zero(dim, numUncertainVariables_);
        tempScalar = alpha_ * std::sqrt(numUncertainVariables_ + kappa_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = tempScalar;
        }
        normalizedPoints_ << Vector::Zero(dim), A, -A;

        // Set weights
        ctypeRNum weight = 1.0 / (2.0 * alpha_*alpha_ * (numUncertainVariables_ + kappa_));
        weightsMean_ << 1.0 - numUncertainVariables_ / (alpha_*alpha_ * (numUncertainVariables_ + kappa_)), Vector::Constant(2*numUncertainVariables_, weight);
        weightsVar_ << weightsMean_(0) + 1.0 - alpha_*alpha_ + beta_, Vector::Constant(2*numUncertainVariables_, weight);
    }

    UnscentedTransformation::UnscentedTransformation(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa)
    : UnscentedTransformation(dim, alpha, beta, kappa, std::vector<bool>(dim, true))
    {
    }

    const Matrix& UnscentedTransformation::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }
        return points_;
    }

    const Vector& UnscentedTransformation::mean(MatrixConstRef points)
    {
        mean_.noalias() = points * weightsMean_;
        return mean_;
    }

    typeRNum UnscentedTransformation::mean1D(RowVectorConstRef points)
    {
        return points * weightsMean_;
    }

    const Matrix& UnscentedTransformation::covariance(MatrixConstRef points)
    {
        // Compute the mean
        mean_.noalias() = points * weightsMean_;

        // Difference between points and mean
        diffToMean_ = points.colwise() - mean_;

        // Compute and return the covariance matrix
        covariance_.noalias() = diffToMean_ * weightsVar_.asDiagonal() * diffToMean_.transpose();
        return covariance_;
    }

    typeRNum UnscentedTransformation::variance(RowVectorConstRef points)
    {
        tempScalar = points * weightsMean_;
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        return diffToMean1D_ * diffToMean1D_.transpose().cwiseProduct(weightsVar_);
    }

    const Vector& UnscentedTransformation::dmean_dpoints_vec(VectorConstRef vec)
    {
        // out = (w_1*v_1 w_1*v_2  w_1*v_3  ...  w_2*v_1  w_2*v_2   w_2*v_3  ...)^T

        // Set and return output
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dmean_dpoints_vec_.segment(i * dim_, dim_) = weightsMean_(i) * vec;
        }
        return dmean_dpoints_vec_;
    }

    const Vector& UnscentedTransformation::dmean1D_dpoints()
    {
        return weightsMean_;
    }

    const Vector& UnscentedTransformation::dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec)
    {
        /*
        * Compute out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        * with dcov_ij/dx = sum_k(weightsVar_k * (x_ik - m_i) * (dx_jk/dx - dm_j/d_x) +
        *                         weightsVar_k * (x_jk - m_j) * (dx_ik/dx - dm_i/d_x))
        */

        // Mapping of the vec input to a matrix corresponding to the covariance matrix
        Eigen::Map<const Matrix> vecMatrix(vec.data(), dim_, dim_);

        // Compute the mean and the difference from each sigma point to the mean
        mean_.noalias() = points * weightsMean_;
        diffToMean_ = points.colwise() - mean_;

        // out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        // Initialize output
        dcov_dpoints_vec_.setZero();             
        
        // Compute dcov_ij/dx
        for(typeInt j = 0; j < dim_; ++j)
        {
            // diffToMean_jk_dx = "1 in element jk" + "weightsMean in row j"
            diffToMean_jk_dx_.row(j) = weightsMean_.transpose();

            for(typeInt i = 0; i < dim_; ++i)
            {
                // diffToMean_ik_dx = "1 in element ik" + "weightsMean in row i"
                diffToMean_ik_dx_.row(i) = weightsMean_.transpose();
                
                // Set dcov to zero
                dcov_.setZero();

                // Go through all sigma points
                for(typeInt k = 0; k < numPoints_; ++k)
                {
                    // Add 1 at elements jk and ik
                    diffToMean_jk_dx_(j, k) += 1.0;
                    diffToMean_ik_dx_(i, k) += 1.0;

                    // Add dcov_ij_dx for sigmapoint k
                    dcov_ += weightsVar_(k) * (diffToMean_(i, k) * diffToMean_jk_dx_ + diffToMean_(j, k) * diffToMean_ik_dx_);

                    // Remove 1 at elements jk and ik
                    diffToMean_jk_dx_(j, k) -= 1.0;
                    diffToMean_ik_dx_(i, k) -= 1.0;
                }
                // Remove row j
                diffToMean_jk_dx_.row(i).setZero();

                // Add vecMatrix_ij * dcov_ij/dx to the output
                dcov_dpoints_vec_ += vecMatrix(i, j) * dcov_.reshaped(dim_ * numPoints_, 1);
            }
            // Remove row i
            diffToMean_ik_dx_.row(j).setZero();
        }
        return dcov_dpoints_vec_;
    }

    const Vector& UnscentedTransformation::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points * weightsMean_;

        // Difference of each point to the mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        
        // derivative of the variance
        tempScalar = diffToMean1D_ * weightsVar_;   
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dvar_dpoints_(i) = 2.0 * (weightsVar_(i) * diffToMean1D_(i) - tempScalar * weightsMean_(i));
        }
        return dvar_dpoints_;
    }

    typeInt UnscentedTransformation::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr UT(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new UnscentedTransformation(dim, alpha, beta, kappa, considerUncertain));
    }

    PointTransformationPtr UT(typeInt dim, typeRNum alpha, typeRNum beta, typeRNum kappa)
    {
        return PointTransformationPtr(new UnscentedTransformation(dim, alpha, beta, kappa));
    }
}