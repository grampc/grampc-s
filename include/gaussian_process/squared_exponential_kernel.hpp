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

#ifndef SQUARED_EXPONENTIAL_KERNEL
#define SQUARED_EXPONENTIAL_KERNEL

#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    class SquaredExponentialKernel : public StationaryKernel
    {
    public:
        // Constructor of a squared exponential kernel k(tau) = sigma^2 * exp(-0.5 * sum_i(tau_i^2 / l_i^2)) with tau = x1 - x2, output standard deviation sigma, and a vector of length-scales l
        SquaredExponentialKernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale);

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const override;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const override;

        // return the input dimension
        virtual typeInt inputDimension() const override;

    private:
        typeInt dimInput_;
        typeRNum sigmaSquared_;
        Vector lengthScaleSquared_;
    };

    // Constructor of a squared exponential kernel k(tau) = sigma^2 * exp(-0.5 * sum_i(tau_i^2 / l_i^2)) with tau = x1 - x2, output standard deviation sigma, and a vector of length-scales l
    StationaryKernelPtr SquaredExponential(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale);
}

#endif // SQUARED_EXPONENTIAL_KERNEL