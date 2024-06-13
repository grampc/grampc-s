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


#ifndef LOCALLY_PERIODIC_KERNEL
#define LOCALLY_PERIODIC_KERNEL

#include "gaussian_process/stationary_kernel.hpp"

namespace grampc
{
    class LocallyPeriodicKernel : public StationaryKernel
    {
    public:
        // Constructor of a locally periodic kernel k(tau) = sigma^2 * exp(-2 * sum_i(sin^2(pi * tau_i / p_i) / l_i^2)) * exp(-0.5 * sum_i(tau_i^2 / l_i^2)) with tau = x1 - x2, output standard deviation sigma, a vector of length-scales l, and a vector of periods p.
        LocallyPeriodicKernel(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale, const std::vector<typeRNum>& period);

        // evaluate the kernel at tau = x1 - x2
        virtual typeRNum evaluate(VectorConstRef tau) const override;

        // gradient of the kernel with respect to the specified elements of tau at the evaluation point tauEval
        virtual void gradient(VectorRef out, VectorConstRef tauEval, IntVectorConstRef derivativeIndices) const override;

        // return the input dimension
        virtual typeInt inputDimension() const override;

    private:
        typeInt dimInput_;
        typeRNum sigmaSquared_;
        Vector lengthScaleSquared_;
        Vector period_;
    };

    // Constructor of a locally periodic kernel k(tau) = sigma^2 * exp(-2 * sum_i(sin^2(pi * tau_i / p_i) / l_i^2)) * exp(-0.5 * sum_i(tau_i^2 / l_i^2)) with tau = x1 - x2, output standard deviation sigma, a vector of length-scales l, and a vector of periods p.
    StationaryKernelPtr LocallyPeriodic(typeInt dimInput, typeRNum sigma, const std::vector<typeRNum>& lengthScale, const std::vector<typeRNum>& period);
}

#endif // LOCALLY_PERIODIC_KERNEL