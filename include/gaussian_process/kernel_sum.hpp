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


#ifndef KERNEL_SUM
#define KERNEL_SUM

#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    class KernelSum : public StationaryKernel
    {
    public:
        // Constructor of a sum of kernels with the same input variables
        KernelSum(const std::vector<StationaryKernelPtr>& kernels);

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const override;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const override;

        // return the input dimension
        virtual typeInt inputDimension() const override;

    private:
        typeInt numKernels_;
        typeInt dimInput_;
        std::vector<StationaryKernelPtr> kernels_;
        mutable Vector temp_vec_dimInput_;
    };

    // Constructor of a sum of kernels with the same input variables
    StationaryKernelPtr SumOfKernels(const std::vector<StationaryKernelPtr>& kernels);
}

#endif // KERNEL_SUM