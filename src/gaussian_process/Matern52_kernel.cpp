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


#include "gaussian_process/Matern52_kernel.hpp"

namespace grampc
{
    Matern52Kernel::Matern52Kernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    : dimInput_(dimInput),
      sigmaSquared_(sigma * sigma),
      lengthScaleSquared_(dimInput)
    {
        for(typeInt i = 0; i < dimInput; ++i)
        {
            lengthScaleSquared_(i) = lengthScale[i] * lengthScale[i];
        }
    }

    typeRNum Matern52Kernel::evaluate(VectorConstRef tau) const
    {
        // compute Euclidean distance
        typeRNum d;
        typeRNum d2 = 0.0;
        for(typeInt i = 0; i < dimInput_; ++i)
        {
           d2 += tau(i) * tau(i) / lengthScaleSquared_(i);
        }
        d = std::sqrt(d2);
        
        return sigmaSquared_ * (1.0 + std::sqrt(5.0) * d + 5.0/3.0 * d2) * std::exp(- std::sqrt(5.0) * d);
    }

    void Matern52Kernel::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        // compute Euclidean distance
        typeRNum d = 0.0;
        for(typeInt i = 0; i < dimInput_; ++i)
        {
           d += tauEval(i) * tauEval(i) / lengthScaleSquared_(i);
        }
        d = std::sqrt(d);

        // compute gradient
        typeRNum temp = - 5.0/3.0 * sigmaSquared_ * std::exp(- std::sqrt(5.0) * d) * (1 + std::sqrt(5.0) * d);
        for(typeInt i = 0; i < derivativeIndices.rows(); ++i)
        {
            out(i) = temp * tauEval(derivativeIndices(i)) / lengthScaleSquared_(derivativeIndices(i));
        }
    }

    typeInt Matern52Kernel::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr Matern52(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    {
        return StationaryKernelPtr(new Matern52Kernel(dimInput, sigma, lengthScale));
    }
}
