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

#ifndef STATIONARY_KERNEL_HPP
#define STATIONARY_KERNEL_HPP

#include "util/grampc_s_constants.hpp"

namespace grampc
{
    // Stationary covariance function of a Gaussian process
    class StationaryKernel
    {
    public:
        virtual ~StationaryKernel() {}

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const = 0;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const = 0;

        // return the input dimension
        virtual typeInt inputDimension() const = 0;
    };

    // Alias
    typedef std::shared_ptr<StationaryKernel> StationaryKernelPtr;
    typedef std::shared_ptr<const StationaryKernel> StationaryKernelConstPtr;
}

#endif // STATIONARY_KERNEL_HPP
