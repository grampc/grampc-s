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

#include "double_integrator_GP_problem_description.hpp"

DoubleIntegratorGPProblemDescription::DoubleIntegratorGPProblemDescription(const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
    : pCost_(pCost), pCon_(pCon)
{
}

void DoubleIntegratorGPProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 2;
	*Nu = 1;
	*Np = 0;
	*Nh = 1;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void DoubleIntegratorGPProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0;
	out[1] = 0;
}

void DoubleIntegratorGPProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
	out[1] = 0;
}

void DoubleIntegratorGPProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
}

void DoubleIntegratorGPProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}

void DoubleIntegratorGPProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    ctypeRNum *udes = param->udes;
    out[0] = pCost_[0] * (u[0] - udes[0]) * (u[0] - udes[0]) + pCost_[1] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[2] * (x[1] - xdes[1]) * (x[1] - xdes[1]);
}

void DoubleIntegratorGPProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = 2 * pCost_[1] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[2] * (x[1] - xdes[1]);
}

void DoubleIntegratorGPProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *udes = param->udes;
	out[0] = 2 * pCost_[0] * (u[0] - udes[0]);
}

void DoubleIntegratorGPProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = pCost_[3] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[4] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +  pCost_[5] * t;
}

void DoubleIntegratorGPProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = 2 * pCost_[3] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[4] * (x[1] - xdes[1]);
}

void DoubleIntegratorGPProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = pCost_[5];
}

void DoubleIntegratorGPProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = -x[1] + pCon_[0]; 
}

void DoubleIntegratorGPProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
    out[1] = -vec[0];
}

void DoubleIntegratorGPProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
}

void DoubleIntegratorGPProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}