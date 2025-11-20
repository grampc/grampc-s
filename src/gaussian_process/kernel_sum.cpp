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


#include "gaussian_process/kernel_sum.hpp"

namespace grampc
{
    KernelSum::KernelSum(const std::vector<StationaryKernelPtr>& kernels)
    : dimInput_(kernels[0]->inputDimension()),
      kernels_(kernels),
      numKernels_((typeInt) kernels.size()),
      temp_vec_dimInput_(dimInput_)
    {
        // Input variables of all kernels must be the same!
    }

    typeRNum KernelSum::evaluate(VectorConstRef tau) const
    {
        typeRNum temp = 0.0;
        for(typeInt i = 0; i < numKernels_; ++i)
        {
            temp += kernels_[i]->evaluate(tau);
        }
        
        return temp;
    }

    void KernelSum::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        out.setZero();
        typeInt numDerivatives = (typeInt) derivativeIndices.rows();

        for (typeInt i = 0; i < numKernels_; ++i)
        {
            kernels_[i]->gradient(temp_vec_dimInput_, tauEval, derivativeIndices);
            out += temp_vec_dimInput_.topRows(numDerivatives);
        }
    }

    typeInt KernelSum::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr SumOfKernels(const std::vector<StationaryKernelPtr>& kernels)
    {
        return StationaryKernelPtr(new KernelSum(kernels));
    }
}
