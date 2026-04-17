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

#include "double_integrator_problem_description.hpp"

DoubleIntegratorProblemDescription::DoubleIntegratorProblemDescription(const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
    : ProblemDescription(2, 1, 0, 0, 1, 0, 0), pCost_(pCost), pCon_(pCon)
{
}

void DoubleIntegratorProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    out[0] = x[1];
	out[1] = u[0];
}

void DoubleIntegratorProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = 0;
	out[1] = vec[0];
}

void DoubleIntegratorProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = vec[1];
}

void DoubleIntegratorProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
}

void DoubleIntegratorProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    auto& xdes = param.xdes;
    auto& udes = param.udes;
    out[0] = pCost_[0] * (u[0] - udes[0]) * (u[0] - udes[0]) + pCost_[1] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[2] * (x[1] - xdes[1]) * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    auto& xdes = param.xdes;
    out[0] = 2 * pCost_[1] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[2] * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    auto& udes = param.udes;
	out[0] = 2 * pCost_[0] * (u[0] - udes[0]);
}

void DoubleIntegratorProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{
    auto& xdes = param.xdes;
    out[0] = pCost_[3] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[4] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +  pCost_[5] * t;
}

void DoubleIntegratorProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{
    auto& xdes = param.xdes;
    out[0] = 2 * pCost_[3] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[4] * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const GrampcParam& param)
{
    out[0] = pCost_[5];
}

void DoubleIntegratorProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const GrampcParam& param)
{
    out[0] = -x[1] + pCon_[0]; 
}

void DoubleIntegratorProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = 0;
    out[1] = -vec[0];
}

void DoubleIntegratorProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
    out[0] = 0;
}

void DoubleIntegratorProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const GrampcParam& param)
{
}