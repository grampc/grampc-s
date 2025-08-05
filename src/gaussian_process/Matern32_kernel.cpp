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


#include "gaussian_process/Matern32_kernel.hpp"

namespace grampc
{
    Matern32Kernel::Matern32Kernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    : dimInput_(dimInput),
      sigmaSquared_(sigma * sigma),
      lengthScaleSquared_(dimInput)
    {
        for(typeInt i = 0; i < dimInput; ++i)
        {
            lengthScaleSquared_(i) = lengthScale[i] * lengthScale[i];
        }
    }

    typeRNum Matern32Kernel::evaluate(VectorConstRef tau) const
    {
        // compute Euclidean distance
        typeRNum d = 0.0;
        for(typeInt i = 0; i < dimInput_; ++i)
        {
           d += tau(i) * tau(i) / lengthScaleSquared_(i);
        }
        d = std::sqrt(d);
        
        return sigmaSquared_ * (1.0 + std::sqrt(3.0) * d) * std::exp(- std::sqrt(3.0) * d);
    }

    void Matern32Kernel::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        // compute Euclidean distance
        typeRNum d = 0.0;
        for(typeInt i = 0; i < dimInput_; ++i)
        {
           d += tauEval(i) * tauEval(i) / lengthScaleSquared_(i);
        }
        d = std::sqrt(d);

        // compute gradient
        typeRNum temp = - 3.0 * sigmaSquared_ * std::exp(- std::sqrt(3.0) * d);
        for(typeInt i = 0; i < derivativeIndices.rows(); ++i)
        {
            out(i) = temp * tauEval(derivativeIndices(i)) / lengthScaleSquared_(derivativeIndices(i));
        }
    }

    typeInt Matern32Kernel::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr Matern32(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    {
        return StationaryKernelPtr(new Matern32Kernel(dimInput, sigma, lengthScale));
    }
}
