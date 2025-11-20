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


#include "gaussian_process/kernel_product.hpp"

namespace grampc
{
    KernelProduct::KernelProduct(const StationaryKernelPtr& kernel1, const StationaryKernelPtr& kernel2)
    : dimInput_(kernel1->inputDimension()),
      kernel1_(kernel1),
      kernel2_(kernel2),      
      temp_vec_dimInput_(dimInput_)
    {
        // Input variables of all kernels must be the same!
    }

    typeRNum KernelProduct::evaluate(VectorConstRef tau) const
    {
        return kernel1_->evaluate(tau) * kernel2_->evaluate(tau);
    }

    void KernelProduct::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        typeInt numDerivatives = (typeInt) derivativeIndices.rows();

        // dk/dx = dk1/dx * k2 + k1 * dk2/dx 
        kernel1_->gradient(temp_vec_dimInput_, tauEval, derivativeIndices);
        out = temp_vec_dimInput_.topRows(numDerivatives) * kernel1_->evaluate(tauEval);
        kernel2_->gradient(temp_vec_dimInput_, tauEval, derivativeIndices);
        out += temp_vec_dimInput_.topRows(numDerivatives) * kernel2_->evaluate(tauEval);
    }

    typeInt KernelProduct::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr ProductOfKernels(const StationaryKernelPtr& kernel1, const StationaryKernelPtr& kernel2)
    {
        return StationaryKernelPtr(new KernelProduct(kernel1, kernel2));
    }
}
