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


#ifndef MATERN_32_KERNEL
#define MATERN_32_KERNEL

#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    class Matern32Kernel : public StationaryKernel
    {
    public:
        // Constructor of a Matern-3/2 kernel k(tau) = sigma^2 * (1 + sqrt(3) * d) * exp(-sqrt(3) * d) with d = sqrt(sum_i(tau_i^2 / l_i^2)), tau = x1 - x2, output standard deviation sigma, a vector of length-scales l.
        Matern32Kernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale);

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const override;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const override;

        // return the input dimension
        virtual typeInt inputDimension() const override;

    private:
        const typeInt dimInput_;
        typeRNum sigmaSquared_;
        Vector lengthScaleSquared_;
    };

    // Constructor of a Matern-3/2 kernel k(tau) = sigma^2 * (1 + sqrt(3) * d) * exp(-sqrt(3) * d) with d = sqrt(sum_i(tau_i^2 / l_i^2)), tau = x1 - x2, output standard deviation sigma, a vector of length-scales l.
    StationaryKernelPtr Matern32(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale);
}

#endif // MATERN_32_KERNEL