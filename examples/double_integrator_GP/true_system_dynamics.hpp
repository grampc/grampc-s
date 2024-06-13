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

#ifndef TRUE_SYSTEM_DYNAMICS_DOUBLE_INTEGRATOR_HPP
#define TRUE_SYSTEM_DYNAMICS_DOUBLE_INTEGRATOR_HPP

#include "util/grampc_s_constants.hpp"

using namespace grampc;

void doubleIntegratorDynamics(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p);

#endif // TRUE_SYSTEM_DYNAMICS_DOUBLE_INTEGRATOR_HPP