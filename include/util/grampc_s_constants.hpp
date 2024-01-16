#ifndef GRAMPC_S_CONSTANTS_HPP
#define GRAMPC_S_CONSTANTS_HPP

#include <random>
#include <Eigen/Dense>
#include "problem_description/grampc_interface.hpp"


namespace grampc
{
    // mathematical constants
    const typeRNum EULER_MASCHERONI_CONSTANT = 0.57721566490153286060;
    const typeRNum PI = 3.14159265358979323846;

    // minimum standard deviation of constraints used to calculate the gradients in dhdx_vec
    const typeRNum stdDevMin_ = 1e-12;

    typedef std::mt19937 RandomNumberGenerator;

    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<typeInt, Eigen::Dynamic, 1> IntVector;
    typedef Eigen::Matrix<typeRNum, 1, Eigen::Dynamic> RowVector;
    typedef Eigen::DiagonalMatrix<typeRNum, Eigen::Dynamic> DiagonalMatrix;
    
    typedef Eigen::Ref<Matrix> MatrixRef;
    typedef Eigen::Ref<Vector> VectorRef;
    typedef Eigen::Ref<IntVector> IntVectorRef;
    typedef Eigen::Ref<RowVector, 0, Eigen::InnerStride<>> RowVectorRef;
    typedef const Eigen::Ref<const Matrix>& MatrixConstRef;
    typedef const Eigen::Ref<const Vector>& VectorConstRef;
    typedef const Eigen::Ref<const IntVector>& IntVectorConstRef;
    typedef const Eigen::Ref<const RowVector, 0, Eigen::InnerStride<>>& RowVectorConstRef;


    enum class PolynomialFamily
    {
        NONE,
        HERMITE,
        LEGENDRE
    };
}

#endif // GRAMPC_S_CONSTANTS_HPP