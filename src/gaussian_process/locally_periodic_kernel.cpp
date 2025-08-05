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


#include "gaussian_process/locally_periodic_kernel.hpp"

namespace grampc
{
    LocallyPeriodicKernel::LocallyPeriodicKernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale, const std::vector<typeRNum>& period)
    : dimInput_(dimInput),
      sigmaSquared_(sigma * sigma),
      lengthScaleSquared_(dimInput),
      period_(dimInput)
    {
        for(typeInt i = 0; i < dimInput; ++i)
        {
            lengthScaleSquared_(i) = lengthScale[i] * lengthScale[i];
            period_(i) = period[i];
        }
    }

    typeRNum LocallyPeriodicKernel::evaluate(VectorConstRef tau) const
    {
        typeRNum temp1 = 0.0;
        typeRNum temp2 = 0.0;
        typeRNum s;

        for(typeInt i = 0; i < dimInput_; ++i)
        {
            s = std::sin(PI * tau(i) / period_(i));
            temp1 +=  s * s / lengthScaleSquared_(i);
            temp2 += tau(i) * tau(i) / lengthScaleSquared_(i);
        }
        
        return sigmaSquared_ * std::exp(-2.0 * temp1) * std::exp(- 0.5 * temp2);
    }

    void LocallyPeriodicKernel::gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const
    {
        typeRNum arg;
        typeRNum kernel = this->evaluate(tauEval);

        for(typeInt i = 0; i < derivativeIndices.rows(); ++i)
        {
            arg = PI * tauEval(derivativeIndices(i)) / period_(derivativeIndices(i));
            out(i) = -kernel / lengthScaleSquared_(derivativeIndices(i)) * (4.0 * PI / period_(derivativeIndices(i)) * std::sin(arg) * std::cos(arg) + tauEval(derivativeIndices(i)));
        }
    }

    typeInt LocallyPeriodicKernel::inputDimension() const
    {
        return dimInput_;
    }

    StationaryKernelPtr LocallyPeriodic(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale, const std::vector<typeRNum>& period)
    {
        return StationaryKernelPtr(new LocallyPeriodicKernel(dimInput, sigma, lengthScale, period));
    }
}
