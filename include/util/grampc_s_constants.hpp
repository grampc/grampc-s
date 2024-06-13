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


#ifndef GRAMPC_S_CONSTANTS_HPP
#define GRAMPC_S_CONSTANTS_HPP

#include <random>
#include <Eigen/Dense>
#include <memory>

extern "C"
{
#include "grampc.h"
}

namespace grampc
{
    // mathematical constants
    constexpr typeRNum EULER_MASCHERONI_CONSTANT = 0.57721566490153286060;
    constexpr typeRNum PI = 3.14159265358979323846;

    // minimum standard deviation of constraints used to calculate the gradients in dhdx_vec
    constexpr typeRNum stdDevMin_ = 1e-12;

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

    // A system function
    typedef std::function<void (VectorRef, ctypeRNum, VectorConstRef, VectorConstRef, VectorConstRef)> SystemFct;

    enum class PolynomialFamily
    {
        NONE,
        HERMITE,
        LEGENDRE
    };

    struct GaussianProcessData
    {
        typeInt inputDimension;
        typeInt numberOfDataPoints;
        typeRNum outputNoiseVariance;
        Matrix inputData;
        Vector outputData;
    };
    
}

#endif // GRAMPC_S_CONSTANTS_HPP