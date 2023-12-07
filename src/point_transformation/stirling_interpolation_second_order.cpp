#include "point_transformation/stirling_interpolation_second_order.hpp"

namespace grampc
{
    StirlingInterpolationSecondOrder::StirlingInterpolationSecondOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    : dim_(dim),
      stepsize_(stepSize),
      stepSizeSquared_(stepSize*stepSize),
      numUncertainVariables_(std::count(considerUncertain.begin(), considerUncertain.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dim, numPoints_),
      points_(dim, numPoints_),
      weightsMean_(2),
      weightsCov_(2),
      mean_(dim),
      covariance_(dim, dim),
      dmean_dpoints_vec_(dim*numPoints_),
      dmean1D_dpoints_(numPoints_),
      dcov_dpoints_vec_(dim*numPoints_),
      dvar_dpoints_(numPoints_),
      tempVecDim_(dim),
      tempVecDim2_(dim)
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

        // Index for uncertain variables
        typeInt uncertIndex = 0;

        // Fill vector with indices of uncertain variables
        for(typeInt i = 0; i < dim; ++i)
        {
            if(considerUncertain[i])
            {
                uncertainIndices[uncertIndex++] = i;
            }
        }

        // weights for mean and covariance computation
        weightsMean_(0) = (stepSizeSquared_ - numUncertainVariables_)/(stepSizeSquared_);
        weightsMean_(1) = 1.0 / (2.0 * stepSizeSquared_);
        weightsCov_(0) = 1.0 / (4.0 * stepSizeSquared_);
        weightsCov_(1) = (stepSizeSquared_ - 1.0) / (4.0 * stepSizeSquared_ * stepSizeSquared_);

        // Derivative of the one dimensional mean with respect to the points;
        dmean1D_dpoints_(0) = weightsMean_(0);
        dmean1D_dpoints_.segment(1, numPoints_ - 1).setConstant(weightsMean_(1));

        // Set normalized Points
        Matrix A = Matrix::Zero(dim, numUncertainVariables_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = stepsize_;
        }
        normalizedPoints_ << Vector::Zero(dim), A, -A;
    }

     StirlingInterpolationSecondOrder::StirlingInterpolationSecondOrder(typeInt dim, typeRNum stepSize)
     : StirlingInterpolationSecondOrder(dim, stepSize, std::vector<bool>(dim, true))
    {
    }

    const Matrix& StirlingInterpolationSecondOrder::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }
        return points_;
    }

    const Vector& StirlingInterpolationSecondOrder::mean(MatrixConstRef points)
    {
        mean_ = weightsMean_(0) * points.col(0);
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            mean_ += weightsMean_(1) * (points.col(i) + points.col(i + numUncertainVariables_));
        }
        return mean_;
    }

     typeRNum StirlingInterpolationSecondOrder::mean1D(RowVectorConstRef points)
    {
        tempScalar_ = weightsMean_(0) * points(0);
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar_ += weightsMean_(1) * (points(i) + points(i + numUncertainVariables_));
        }
        return tempScalar_;
    }

    const Matrix& StirlingInterpolationSecondOrder::covariance(MatrixConstRef points)
    {
        covariance_.setZero();
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempVecDim_ = points.col(i) - points.col(i + numUncertainVariables_);
            tempVecDim2_ = points.col(i) + points.col(i + numUncertainVariables_) - 2.0 * points.col(0);
            covariance_.noalias() += tempVecDim_ * tempVecDim_.transpose() * weightsCov_(0) + tempVecDim2_ * tempVecDim2_.transpose() * weightsCov_(1);
        }
        return covariance_;
    }

    typeRNum StirlingInterpolationSecondOrder::variance(RowVectorConstRef points)
    {
        tempScalar_ = 0.0;
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar2_ = points(i) - points(i + numUncertainVariables_);
            tempScalar3_ = points(i) + points(i + numUncertainVariables_) - 2.0 * points(0);
            tempScalar_ += tempScalar2_ * tempScalar2_ * weightsCov_(0) + tempScalar3_ * tempScalar3_ * weightsCov_(1);
        }
        return tempScalar_;
    }

     const Vector& StirlingInterpolationSecondOrder::dmean_dpoints_vec(VectorConstRef vec)
    {
        dmean_dpoints_vec_.segment(0, dim_) = weightsMean_(0) * vec;

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            dmean_dpoints_vec_.segment(i*dim_, dim_) = weightsMean_(1) * vec;
            dmean_dpoints_vec_.segment((i + numUncertainVariables_) * dim_, dim_) = weightsMean_(1) * vec;
        }
        return dmean_dpoints_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dmean1D_dpoints()
    {
        return dmean1D_dpoints_;
    }

    const Vector& StirlingInterpolationSecondOrder::dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec)
    {
        std::cerr << "not implemented!" << std::endl;
        return dcov_dpoints_vec_;
    }

    const Vector& StirlingInterpolationSecondOrder::dvar_dpoints(RowVectorConstRef points)
    {
        tempScalar_ = 1.0/(2.0 * stepsize_ * stepsize_);
        tempScalar2_ = 2.0 * weightsCov_(1);

        dvar_dpoints_(0) = 0.0;

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar3_ = tempScalar2_ * (points(i) + points(i + numUncertainVariables_) - 2.0 * points(0));
            tempScalar4_ = tempScalar_ * (points(i) - points(i + numUncertainVariables_));
            dvar_dpoints_(i) = tempScalar3_ + tempScalar4_;
            dvar_dpoints_(i + numUncertainVariables_) = tempScalar3_ - tempScalar4_;
            dvar_dpoints_(0) -= 2.0 * tempScalar3_;
        }
        return dvar_dpoints_;
    }

     typeInt StirlingInterpolationSecondOrder::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr StirlingSecondOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new StirlingInterpolationSecondOrder(dim, stepSize, considerUncertain));
    }

    PointTransformationPtr StirlingSecondOrder(typeInt dim, typeRNum stepSize)
    {
        return PointTransformationPtr(new StirlingInterpolationSecondOrder(dim, stepSize));
    }
}