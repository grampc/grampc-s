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


#include "point_transformation/composed_quadrature.hpp"

namespace grampc
{
    ComposedQuadrature::ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
      : dimX_(dimX),
        dimY_(dimY),
        numPoints_(quadratureOrder.prod()),
        normalizedPoints_(dimX_, numPoints_),
        points_(dimX_, numPoints_),
        roots_(dimX_, numPoints_),
        meanX_(dimX_),
        meanY_(dimY_),
        covariance_(dimX_, dimY_),
        weights_(Vector::Ones(numPoints_)),
        dvar_dpoints_(numPoints_),
        diffToMeanX_(dimX_, numPoints_),
        diffToMeanY_(dimY_, numPoints_),
        diffToMean1D_(numPoints_)
    {   
        // Check inputs 
        if(polyFamily.size() != quadratureOrder.size())
        {
            std::cerr << "Dimensions of inputs do not match!" << std::endl;
        }

        // Get corresponding quadrature rules
        std::vector<QuadratureRuleConstPtr> quadRules(dimX_);
        for(typeInt i = 0; i < dimX_; ++i)
        {
            quadRules[i] = correspondingQuadratureRule(polyFamily[i], quadratureOrder[i]);
        }

        /************************************************************************************
         *  Initialize roots, normalized points, and weights
         ************************************************************************************/
        typeInt index;

        // Vector of quadrature points numbers for the univariate quadrature rules
        std::vector<typeInt> numUnivariateQuadPoints(dimX);

        // Number of quadrature points for the multivariate quadrature
        typeInt numOutputPoints = 1;

        // Compute vector of quadrature points numbers and number of output points
        for(typeInt i = 0; i < dimX; ++i)
        {
            numUnivariateQuadPoints[i] = quadRules[i]->numberOfPoints();
            numOutputPoints *= numUnivariateQuadPoints[i] ;
        }

        // Number of possible combinations
        typeInt numComb = numOutputPoints;

        // Compute roots and normalized quadrature points
        for(typeInt i = 0; i < dimX; ++i)
        {
            // Get points and weights of the current quadrature rule
            const Vector& roots = quadRules[i]->polynomialRoots();
            const Vector& normalizedPoints = quadRules[i]->pointsNormalized();
            const Vector& weights = quadRules[i]->weights();

            // Add points and update weights
            for(typeInt j = 0; j < numOutputPoints; ++j)
            {
                index = (typeInt) std::floor((j % numComb) / (numComb / numUnivariateQuadPoints[i]));
                roots_(i, j) = roots(index);
                normalizedPoints_(i, j) = normalizedPoints(index);
                weights_(j) *= weights(index);
            }
            numComb /= numUnivariateQuadPoints[i];
        }
    }

    ComposedQuadrature::ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder)
    : ComposedQuadrature(dimX, dimY, polyFamily, Eigen::Map<const Eigen::Vector<typeInt, Eigen::Dynamic>>(quadratureOrder.data(), polyFamily.size()))
    {
    }

    ComposedQuadrature::ComposedQuadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder)
    : ComposedQuadrature(dimX, dimY, polyFamily, Eigen::Vector<typeInt, Eigen::Dynamic>::Constant(polyFamily.size(), quadratureOrder))
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

    const Matrix& ComposedQuadrature::points(VectorConstRef mean, MatrixConstRef covCholesky)
    {
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = mean;
            points_.col(i).noalias() += covCholesky * normalizedPoints_.col(i);
        }

        return points_;
    }

    const Matrix& ComposedQuadrature::points()
    {
        return points_;
    }

    const Vector& ComposedQuadrature::mean(MatrixConstRef points)
    {
        meanY_.noalias() = points * weights_;
        return meanY_;
    }

    typeRNum ComposedQuadrature::mean1D(RowVectorConstRef points)
    {
        return points * weights_;
    }

    const Matrix& ComposedQuadrature::covariance(MatrixConstRef pointsX, MatrixConstRef pointsY)
    {
        // Compute mean of points
        meanX_.noalias() = pointsX * weights_;
        meanY_.noalias() = pointsY * weights_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMeanX_.col(i) = (pointsX.col(i) - meanX_) * weights_(i);
            diffToMeanY_.col(i) = pointsY.col(i) - meanY_;
        }

        // Compute and return the covariance matrix
        covariance_.noalias() = diffToMeanX_ * diffToMeanY_.transpose();
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

    const Vector& ComposedQuadrature::dmean1D_dpoints()
    {
        return weights_;
    }

    const Vector& ComposedQuadrature::dvar_dpoints(RowVectorConstRef points)
    {
        // Compute mean
        tempScalar = points * weights_;

        // difference between points and mean
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            diffToMean1D_(i) = points(i) - tempScalar;
        }

        tempScalar = - 2.0 * diffToMean1D_ * weights_.cwiseProduct(weights_); 
        
        // Derivative of the variance
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            dvar_dpoints_(i) = 2.0 * weights_(i) * diffToMean1D_(i) + tempScalar;
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

    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>& quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(dimX, dimY, polyFamily, quadratureOrder));
    }

    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, const std::vector<typeInt>& quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(dimX, dimY, polyFamily, quadratureOrder));
    }

    PointTransformationPtr Quadrature(typeInt dimX, typeInt dimY, const std::vector<PolynomialFamily>& polyFamily, typeInt quadratureOrder)
    {
        return PointTransformationPtr(new ComposedQuadrature(dimX, dimY, polyFamily, quadratureOrder));
    }

} // namespace grampc
