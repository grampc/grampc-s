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


#include "gaussian_process/squared_exponential_kernel.hpp"

namespace grampc
{
    SquaredExponentialKernel::SquaredExponentialKernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    : dimInput_(dimInput),
      sigmaSquared_(sigma * sigma),
      lengthScaleSquared_(lengthScale.size())
    {
        for(typeInt i = 0; i < lengthScale.size(); ++i)
        {
            lengthScaleSquared_(i) = lengthScale[i] * lengthScale[i];
        }
    }

    typeRNum SquaredExponentialKernel::evaluate(VectorConstRef tau) const
    {
        typeRNum temp = 0.0;

        for(typeInt i = 0; i < dimInput_; ++i)
        {
            temp += tau(i) * tau(i) / lengthScaleSquared_(i);
        }
        
        return sigmaSquared_ * std::exp(- 0.5 * temp);
    }

    void SquaredExponentialKernel::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        typeRNum kernel = this->evaluate(tauEval);

        for(typeInt i = 0; i < derivativeIndices.rows(); ++i)
        {
            out(i) = -kernel * tauEval(derivativeIndices(i)) / lengthScaleSquared_(derivativeIndices(i));
        }
    }

    typeInt SquaredExponentialKernel::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr SquaredExponential(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale)
    {
        return StationaryKernelPtr(new SquaredExponentialKernel(dimInput, sigma, lengthScale));
    }
}
