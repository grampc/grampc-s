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


#ifndef KERNEL_PRODUCT
#define KERNEL_PRODUCT

#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    class KernelProduct : public StationaryKernel
    {
    public:
        // Constructor of a product of two kernels k1 and k2 with the same input variables
        KernelProduct(const StationaryKernelPtr& kernel1, const StationaryKernelPtr& kernel2);

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const override;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const override;

        // return the input dimension
        virtual typeInt inputDimension() const override;

    private:
        typeInt dimInput_;
        StationaryKernelPtr kernel1_;
        StationaryKernelPtr kernel2_;
        mutable Vector temp_vec_dimInput_;
    };

    // Constructor of a product of two kernels k1 and k2 with the same input variables
    StationaryKernelPtr ProductOfKernels(const StationaryKernelPtr& kernel1, const StationaryKernelPtr& kernel2);
}

#endif // KERNEL_PRODUCT