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


#include "point_transformation/PCE_transformation.hpp"

namespace grampc
{
    PCE_Transformation::PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
    : dimX_(dimX),
      dimY_(dimY),
      dimUncertain_((quadratureOrder.array() > 1).count()),
      numPoints_(quadratureOrder.prod()),
      maxPolyOrder_(maxPolyOrder),
      numPolynomials_(factorial(dimUncertain_ + maxPolyOrder_) / factorial(dimUncertain_) / factorial(maxPolyOrder_)),
      normalizedPoints_(dimX_, numPoints_),
      points_(dimX_, numPoints_),
      PCE_coefficients_mean_(dimY_),
      PCE_coefficientsX_(dimX_, numPolynomials_ - 1),
      PCE_coefficientsY_(dimY_, numPolynomials_ - 1),
      PCE_coefficients_var_(numPolynomials_-1),
      squaredNorms_(numPolynomials_),
      squaredNorms_Var_(numPolynomials_-1),
      temp_squaredNorms_Var_(numPolynomials_-1),
      covariance_(dimX_, dimY_),
      weightsCov_(numPoints_, numPolynomials_-1),
      weightsMean_(numPoints_),
      dvar_dpoints_vec_(numPoints_)
    {
        // Quadrature points
        ComposedQuadrature quad = ComposedQuadrature(dimX, dimY, polyFamily, quadratureOrder);
        normalizedPoints_ = quad.normalizedPoints();

        // Construct a map from polynomial families to generators for the corresponding orthogonal polynomials
        std::map<PolynomialFamily, PolynomialGeneratorPtr> polynomialMap;
        for(PolynomialFamily family : polyFamily)
        {
            // Insert an entry if the polynomial family is not contained
            if(polynomialMap.find(family) == polynomialMap.end())
            {
                polynomialMap.insert({family, correspondingPolynomialGenerator(family, maxPolyOrder_)});
            }
        }

        // Compute Multi-indices for the uncertain variables
        Eigen::Matrix<typeInt, -1, -1> multiIndexUncertain = Eigen::Matrix<typeInt, -1, -1>::Zero(numPolynomials_, dimUncertain_);
        for(typeInt i = 1; i < numPolynomials_; ++i)
        {
            if(multiIndexUncertain.row(i-1).sum() < maxPolyOrder_)
            {
                multiIndexUncertain.block(i, 0, 1, dimUncertain_ - 1) = multiIndexUncertain.block(i - 1, 0, 1, dimUncertain_ - 1);
                multiIndexUncertain(i, dimUncertain_-1) = multiIndexUncertain(i-1, dimUncertain_-1) + 1;
            }else{
                for(typeInt j = dimUncertain_ - 1; j >= 0; --j)
                {
                    if(multiIndexUncertain(i-1, j) > 0)
                    {
                        multiIndexUncertain(i, j-1) = multiIndexUncertain(i-1, j-1) + 1;
                        if(j>1)
                        {
                            multiIndexUncertain.block(i, 0, 1, j-1) = multiIndexUncertain.block(i-1, 0, 1, j-1);
                        }
                        break;
                    }
                }
            }
        }

        // Extend multiindex by inserting non-uncertain variables
        Eigen::Matrix<typeInt, -1, -1> multiIndex(numPolynomials_, dimX_);
        if(dimX_ == dimUncertain_)
        {
            multiIndex = std::move(multiIndexUncertain);
        }else{
            typeInt uncertainVarCounter = 0;
            for(typeInt i = 0; i < dimX_; ++i)
            {
                if(quadratureOrder(i) > 1)
                {
                    multiIndex.col(i) = multiIndexUncertain.col(uncertainVarCounter);
                    uncertainVarCounter++;
                }else{
                    multiIndex.col(i).setZero();
                }
            }
        }
        
        // Vector of univariate and multivariate polynomials
        std::vector<PolynomialConstPtr> uniPolyVec(dimX_);
        std::vector<typeRNum> squaredNorms(dimX_);
        std::vector<MultivariatePolynomialPtr> multPolyVec(numPolynomials_);        

        // Construct multivariate polynomials as products of the univariate polynomilas
        for(typeInt i = 0; i < numPolynomials_; ++i)
        {
            for(typeInt j = 0; j < dimX_; ++j)
            {
                uniPolyVec[j] = polynomialMap[polyFamily[j]]->getPolynomial(multiIndex(i, j));
                squaredNorms[j] = polynomialMap[polyFamily[j]]->getSquaredNorm(multiIndex(i, j));
            }
            multPolyVec[i] = MultivarPoly(uniPolyVec, squaredNorms);
            squaredNorms_[i] = multPolyVec[i]->squaredNorm();
        }

        // Squared norms for the computation of the variance
        squaredNorms_Var_ = squaredNorms_.segment(1, numPolynomials_-1);

        // Compute weights for the computation of the polynomial coefficients for the mean
        for(typeInt j = 0; j < numPoints_; ++j)
        {
                weightsMean_(j) = quad.weights()[j] * multPolyVec[0]->evaluate(quad.roots().col(j))  / squaredNorms_[0];
        }

        // Compute weights for the computation of the polynomial coefficients for the covariance matrix
        for(typeInt i = 0; i < numPolynomials_ - 1; ++i)
        {
            for(typeInt j = 0; j < numPoints_; ++j)
            {
                weightsCov_(j, i) = quad.weights()[j] * multPolyVec[i+1]->evaluate(quad.roots().col(j))  / squaredNorms_Var_[i];
            }
        }
    }

    PCE_Transformation::PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const std::vector<typeInt>& quadratureOrder)
    : PCE_Transformation(dimX, dimY, polyFamily, maxPolyOrder, Eigen::Map<const Eigen::Vector<typeInt, Eigen::Dynamic>>(quadratureOrder.data(), polyFamily.size()))
    {
    }

    PCE_Transformation::PCE_Transformation(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, typeInt quadratureOrder)
    : PCE_Transformation(dimX, dimY, polyFamily, maxPolyOrder, Eigen::Vector<typeInt, Eigen::Dynamic>::Constant(polyFamily.size(), quadratureOrder))
    {
    }

    const Matrix& PCE_Transformation::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();

        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& PCE_Transformation::points()
    {
        return points_;
    }

    const Vector& PCE_Transformation::mean(MatrixConstRef points)
    {
        PCE_coefficients_mean_.noalias() = points * weightsMean_;
        return PCE_coefficients_mean_;
    }

    typeRNum PCE_Transformation::mean1D(RowVectorConstRef points)
    {
        return points * weightsMean_;
    }

    const Matrix& PCE_Transformation::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        // Compute the polynomial coefficients
        PCE_coefficientsX_.noalias() = pointsX * weightsCov_;
        PCE_coefficientsY_.noalias() = pointsY * weightsCov_;

        // cross-covariance between x and y
        for(typeInt i = 0; i < dimX_; ++i)
        {
            for(typeInt j = 0; j < dimY_; ++j)
            {
                covariance_(i, j) = (PCE_coefficientsX_.row(i).cwiseProduct(PCE_coefficientsY_.row(j))  * squaredNorms_Var_).value();
            }
        }
        return covariance_;
    }

    typeRNum PCE_Transformation::variance(RowVectorConstRef points)
    {
        // Compute the polynomial coefficients
        PCE_coefficients_var_.noalias() = points * weightsCov_;
        PCE_coefficients_var_ = PCE_coefficients_var_.cwiseAbs2();
        
        return PCE_coefficients_var_ * squaredNorms_Var_;
    }

    const Vector& PCE_Transformation::dmean1D_dpoints()
    {
        return weightsMean_;
    }

    const Vector& PCE_Transformation::dvar_dpoints(RowVectorConstRef points)
    {
        // Polynomial coefficients
        PCE_coefficients_var_.noalias() = points * weightsCov_;

        // Derivative of the variance
        temp_squaredNorms_Var_ = squaredNorms_Var_.cwiseProduct(PCE_coefficients_var_.transpose());
        dvar_dpoints_vec_.noalias() = 2.0 * weightsCov_ * temp_squaredNorms_Var_;

        return dvar_dpoints_vec_;
    }

    typeInt PCE_Transformation::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
    {
        return PointTransformationPtr(new PCE_Transformation(dimX, dimY, polyFamily, maxPolyOrder, quadratureOrder));
    }

    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, const std::vector<typeInt>& quadratureOrder)
    {
        return PointTransformationPtr(new PCE_Transformation(dimX, dimY, polyFamily, maxPolyOrder, quadratureOrder));
    }

    PointTransformationPtr PCE(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt maxPolyOrder, typeInt quadratureOrder)
    {
        return PointTransformationPtr(new PCE_Transformation(dimX, dimY, polyFamily, maxPolyOrder, quadratureOrder));
    }
}