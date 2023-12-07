#include "point_transformation/Stirling_interpolation_first_order.hpp"


namespace grampc
{
    StirlingInterpolationFirstOrder::StirlingInterpolationFirstOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    : dim_(dim),
      stepsize_(stepSize),
      numUncertainVariables_(std::count(considerUncertain.begin(), considerUncertain.end(), true)),
      numPoints_(2*numUncertainVariables_ + 1),
      normalizedPoints_(dim, numPoints_),
      points_(dim, numPoints_),
      mean_(dim),
      covariance_(dim, dim),
      dmean_dpoints_vec_(dim*numPoints_),
      dcov_dpoints_vec_(dim*numPoints_),
      dvar_dpoints_(numPoints_),
      tempVecDim_(dim)
    {
        // Vector of indices with uncertain variables
        Eigen::Vector<typeInt, Eigen::Dynamic> uncertainIndices(numUncertainVariables_);

        // Derivative of the one dimensional mean with respect to the points;
        dmean1D_dpoints_ = Vector::Zero(numPoints_);
        dmean1D_dpoints_(0) = 1.0;

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

        // Set normalized Points
        Matrix A = Matrix::Zero(dim, numUncertainVariables_);
        for(typeInt i = 0; i < numUncertainVariables_; ++i)
        {
            A(uncertainIndices(i), i) = stepsize_;
        }
        normalizedPoints_ << Vector::Zero(dim), A, -A;
    }

    StirlingInterpolationFirstOrder::StirlingInterpolationFirstOrder(typeInt dim, typeRNum stepSize)
    : StirlingInterpolationFirstOrder(dim, stepSize, std::vector<bool>(dim, true))
    {
    }

    const Matrix& StirlingInterpolationFirstOrder::points(DistributionConstPtr dist)
    {
        const Matrix& chol = dist->covCholesky();
        for(typeInt i = 0; i < numPoints_; ++i)
        {
            points_.col(i) = dist->mean();
            points_.col(i).noalias() += chol * normalizedPoints_.col(i);
        }
        return points_;
    }

    const Vector& StirlingInterpolationFirstOrder::mean(MatrixConstRef points)
    {
        mean_ = points.col(0);
        return mean_;
    }

    typeRNum StirlingInterpolationFirstOrder::mean1D(RowVectorConstRef points)
    {
        return points(0);
    }
    
    const Matrix& StirlingInterpolationFirstOrder::covariance(MatrixConstRef points)
    {
        covariance_.setZero();
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempVecDim_ = points_.col(i) - points.col(i + numUncertainVariables_);
            covariance_.noalias() += tempVecDim_ * tempVecDim_.transpose();
        }
        covariance_ /= (4.0 * stepsize_ * stepsize_);
        return covariance_;
    }

    typeRNum StirlingInterpolationFirstOrder::variance(RowVectorConstRef points)
    {
        tempScalar = 0.0;
        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            tempScalar2 = points(i) - points(i + numUncertainVariables_);
            tempScalar += tempScalar2 * tempScalar2;
        }
        return tempScalar / (4.0 * stepsize_ * stepsize_);
    }

    const Vector& StirlingInterpolationFirstOrder::dmean_dpoints_vec(VectorConstRef vec)
    {
        for(typeInt i = 0; i < dim_; ++i)
        {
            dmean_dpoints_vec_(i) = vec(i);
        }

        for(typeInt i = dim_; i < dim_*numPoints_; ++i)
        {
            dmean_dpoints_vec_(i) = 0.0;
        }

        return dmean_dpoints_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dmean1D_dpoints()
    {
        return dmean1D_dpoints_;
    }

    const Vector& StirlingInterpolationFirstOrder::dcov_dpoints_vec(MatrixConstRef points, VectorConstRef vec)
    {
        std::cerr << "not implemented!" << std::endl;
        return dcov_dpoints_vec_;
    }

    const Vector& StirlingInterpolationFirstOrder::dvar_dpoints(RowVectorConstRef points)
    {
        dvar_dpoints_(0) = 0.0;
        tempScalar = 1.0/(2.0 * stepsize_ * stepsize_);

        for(typeInt i = 1; i <= numUncertainVariables_; ++i)
        {
            dvar_dpoints_(i) = tempScalar * (points(i) - points(i + numUncertainVariables_));
            dvar_dpoints_(i + numUncertainVariables_) = -dvar_dpoints_(i);
        }
        return dvar_dpoints_;
    }

     typeInt StirlingInterpolationFirstOrder::numberOfPoints() const
    {
        return numPoints_;
    }

    PointTransformationPtr StirlingFirstOrder(typeInt dim, typeRNum stepSize, const std::vector<bool>& considerUncertain)
    {
        return PointTransformationPtr(new StirlingInterpolationFirstOrder(dim, stepSize, considerUncertain));
    }

    PointTransformationPtr StirlingFirstOrder(typeInt dim, typeRNum stepSize)
    {
        return PointTransformationPtr(new StirlingInterpolationFirstOrder(dim, stepSize));
    }
}