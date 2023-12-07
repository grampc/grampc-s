#include "point_transformation/composed_quadrature.hpp"

namespace grampc
{
    ComposedQuadrature::ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
      : dim_(polyFamily.size()),
        numPoints_(quadratureOrder.prod()),
        normalizedPoints_(dim_, numPoints_),
        points_(dim_, numPoints_),
        roots_(dim_, numPoints_),
        mean_(dim_),
        covariance_(dim_, dim_),
        covariance_temp_(dim_, dim_),
        weights_(Vector::Ones(numPoints_)),
        dmean_dpoints_vec_(dim_ * numPoints_),
        dcov_dpoints_vec_(dim_ * numPoints_),
        dvar_dpoints_(numPoints_),
        diffToMean_(dim_, numPoints_),
        diffToMean1D_(numPoints_),
        dcov_(dim_, numPoints_),
        diffToMean_jk_dx_(Matrix::Zero(dim_, numPoints_)),
        diffToMean_ik_dx_(Matrix::Zero(dim_, numPoints_))
    {   
        // Check inputs 
        if(polyFamily.size() != quadratureOrder.size())
        {
            std::cerr << "Dimensions of inputs do not match!" << std::endl;
        }

        // Get corresponding quadrature rules
        std::vector<QuadratureRuleConstPtr> quadRules(dim_);
        for(typeInt i = 0; i < dim_; ++i)
        {
            quadRules[i] = correspondingQuadratureRule(polyFamily[i], quadratureOrder[i]);
        }

        /************************************************************************************
         *  Initialize roots, normalized points, and weights
         ************************************************************************************/
        typeInt index;

        // Dimension and vector of quadrature points numbers for the univariate quadrature rules
        typeInt dim = quadRules.size();
        std::vector<typeInt> numUnivariateQuadPoints(dim);

        // Number of quadrature points for the multivariate quadrature
        typeInt numOutputPoints = 1;

        // Compute vector of quadrature points numbers and number of output points
        for(typeInt i = 0; i < dim; ++i)
        {
            numUnivariateQuadPoints[i] = quadRules[i]->numberOfPoints();
            numOutputPoints *= numUnivariateQuadPoints[i] ;
        }

        // Number of possible combinations
        typeInt numComb = numOutputPoints;

        // Compute roots and normalized quadrature points
        for(typeInt i = 0; i < dim; ++i)
        {
            // Get points and weights of the current quadrature rule
            const Vector& roots = quadRules[i]->polynomialRoots();
            const Vector& normalizedPoints = quadRules[i]->pointsNormalized();
            const Vector& weights = quadRules[i]->weights();

            // Add points and update weights
            for(typeInt j = 0; j < numOutputPoints; ++j)
            {
                index = std::floor((j % numComb) / (numComb / numUnivariateQuadPoints[i]));
                roots_(i, j) = roots(index);
                normalizedPoints_(i, j) = normalizedPoints(index);
                weights_(j) *= weights(index);
            }
            numComb /= numUnivariateQuadPoints[i];
        }
    }

    ComposedQuadrature::ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder)
    : ComposedQuadrature(polyFamily, Eigen::Map<const Eigen::Vector<typeInt, Eigen::Dynamic>>(quadratureOrder.data(), polyFamily.size()))
    {
    }

    ComposedQuadrature::ComposedQuadrature(const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder)
    : ComposedQuadrature(polyFamily, Eigen::Vector<typeInt, Eigen::Dynamic>::Constant(polyFamily.size(), quadratureOrder))
    {
    }

    const Matrix& ComposedQuadrature::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }
        return points_;
    }

    const Vector& ComposedQuadrature::mean(MatrixConstRef points)
    {
        mean_.noalias() = points * weights_;
        return mean_;
    }

    typeRNum ComposedQuadrature::mean1D(RowVectorConstRef points)
    {
        return points * weights_;
    }

    const Matrix& ComposedQuadrature::covariance(MatrixConstRef points)
    {
        // Compute the mean
        mean_.noalias() = points * weights_;

        // Difference between points and mean
        diffToMean_ = points.colwise() - mean_;

        // Compute and return the covariance matrix
        covariance_.noalias() = diffToMean_ * weights_.asDiagonal() * diffToMean_.transpose();
        return covariance_;
    }

    typeRNum ComposedQuadrature::variance(RowVectorConstRef points)
    {
        tempScalar = points * weights_;
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        return diffToMean1D_ * diffToMean1D_.transpose().cwiseProduct(weights_);
    }

    const Vector& ComposedQuadrature::dmean_dpoints_vec(VectorConstRef vec)
    {
        // out = (w_1*v_1  w_1*v_2  w_1*v_3  ...  w_2*v_1  w_2*v_2   w_2*v_3  ...)^T

        // Set and return output
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dmean_dpoints_vec_.segment(i * dim_, dim_) = weights_(i) * vec;
        }
        return dmean_dpoints_vec_;
    }

    const Vector& ComposedQuadrature::dmean1D_dpoints()
    {
        return weights_;
    }

    const Vector& ComposedQuadrature::dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec)
    {
        /*
        * Compute out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        * with dcov_ij/dx = sum_k(weightsVar_k * (x_ik - m_i) * (dx_jk/dx - dm_j/d_x) +
        *                         weightsVar_k * (x_jk - m_j) * (dx_ik/dx - dm_i/d_x))
        */

        // Mapping of the vec input to a matrix corresponding to the covariance matrix
        Eigen::Map<const Matrix> vecMatrix(vec.data(), dim_, dim_);

        // Compute the mean and the difference from each sigma point to the mean
        mean_.noalias() = points * weights_;
        diffToMean_ = points.colwise() - mean_;

        // out = vec_1 * dcov_11/dx + vec_2 * dcov_21/dx + ... 
        // Initialize output
        dcov_dpoints_vec_.setZero();             
        
        // Compute dcov_ij/dx
        for(typeInt j = 0; j < dim_; ++j)
        {
            // diffToMean_jk_dx = "1 in element jk" + "weightsMean in row j"
            diffToMean_jk_dx_.row(j) = weights_.transpose();

            for(typeInt i = 0; i < dim_; ++i)
            {
                // diffToMean_ik_dx = "1 in element ik" + "weightsMean in row i"
                diffToMean_ik_dx_.row(i) = weights_.transpose();
                
                // Set dcov to zero
                dcov_.setZero();

                // Go through all sigma points
                for(typeInt k = 0; k < numPoints_; ++k)
                {
                    // Add 1 at elements jk and ik
                    diffToMean_jk_dx_(j, k) += 1.0;
                    diffToMean_ik_dx_(i, k) += 1.0;

                    // Add dcov_ij_dx for sigmapoint k
                    dcov_ += weights_(k) * (diffToMean_(i, k) * diffToMean_jk_dx_ + diffToMean_(j, k) * diffToMean_ik_dx_);

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

    const Vector& ComposedQuadrature::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points * weights_;

        // Difference of each point to the mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }
        
        // derivative of the variance
        tempScalar = diffToMean1D_ * weights_;   
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dvar_dpoints_(i) = 2.0 * weights_(i) * (diffToMean1D_(i) - tempScalar);
        }
        return dvar_dpoints_;
    }

    typeInt ComposedQuadrature::numberOfPoints() const
    {
        return numPoints_;
    }

    const Matrix& ComposedQuadrature::roots()
    {
        return roots_;
    }

    const Vector& ComposedQuadrature::weights()
    {
        return weights_;
    }

    const Matrix& ComposedQuadrature::normalizedPoints()
    {
        return normalizedPoints_;
    }

    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(polyFamily, quadratureOrder));
    }

    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(polyFamily, quadratureOrder));
    }

    PointTransformationPtr Quadrature(const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(polyFamily, quadratureOrder));
    }

} // namespace grampc
