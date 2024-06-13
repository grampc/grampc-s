/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 * 
 * This file describes a modified crane 2D problem from
 * Kapernick, B., Graichen, K.: Model predictive control of an overhead
 * crane using constraint substitution. In: Proc. American Control 
 * Conference (ACC), pp. 3973-3978, 2013.
 * 
 * Based on GRAMPC:
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * 
 */

#include "crane_2d_problem_description.hpp"

 /* square macro */
#define POW2(a) ((a)*(a))

Crane2DProblemDescription::Crane2DProblemDescription(const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
: pCost_(pCost),
  pCon_(pCon)
{
}

void Crane2DProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 6;
	*Nu = 2;
	*Np = 0;
	*Nh = 3;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void Crane2DProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
	out[0] = x[1];
	out[1] = u[0];
	out[2] = x[3];
	out[3] = u[1];
	out[4] = x[5];
	out[5] = -((9.81 * std::sin(x[4]) + std::cos(x[4])*u[0] + 2 * x[3] * x[5]) / x[2]);
}

void Crane2DProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
	ctypeRNum aux1 = std::sin(x[4]);
	ctypeRNum aux2 = std::cos(x[4]);

	out[0] = 0;
	out[1] = vec[0];
	out[2] = (9.81 * aux1 + aux2 * u[0] + 2 * x[3] * x[5])*vec[5] / (x[2] * x[2]);
	out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2];
	out[4] = -(((9.81 * aux2 - aux1 * u[0])*vec[5]) / x[2]);
	out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2];
}

void Crane2DProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
	out[0] = vec[1] - ((std::cos(x[4]))*vec[5]) / x[2];
	out[1] = vec[3];
}

void Crane2DProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
}

void Crane2DProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
	out[0] = pCost_[0] * POW2(x[0] - xdes[0]) +
		pCost_[1] * POW2(x[1] - xdes[1]) +
		pCost_[2] * POW2(x[2] - xdes[2]) +
		pCost_[3] * POW2(x[3] - xdes[3]) +
		pCost_[4] * POW2(x[4] - xdes[4]) +
		pCost_[5] * POW2(x[5] - xdes[5]) +
		pCost_[6] * POW2(u[0] - udes[0]) +
		pCost_[7] * POW2(u[1] - udes[1]);
}

void Crane2DProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
	out[0] = 2 * pCost_[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost_[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost_[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost_[3] * (x[3] - xdes[3]);
	out[4] = 2 * pCost_[4] * (x[4] - xdes[4]);
	out[5] = 2 * pCost_[5] * (x[5] - xdes[5]);
}

void Crane2DProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
	out[0] = 2 * pCost_[6] * (u[0] - udes[0]);
	out[1] = 2 * pCost_[7] * (u[1] - udes[1]);
}

void Crane2DProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
}

void Crane2DProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
}

void Crane2DProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
}

void Crane2DProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
	out[0] = std::cos(x[4]) * x[2] - pCon_[0] * POW2(x[0] + std::sin(x[4])*x[2]) - pCon_[1];
	out[1] = x[5] - pCon_[2];
	out[2] = -x[5] - pCon_[2];
}

void Crane2DProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
	ctypeRNum temp = pCon_[0] * (x[0] + std::sin(x[4])*x[2]);

	out[0] = temp * vec[0];
	out[1] = 0;
	out[2] = (std::sin(x[4])*temp + std::cos(x[4]))*vec[0];
	out[3] = 0;
	out[4] = (std::cos(x[4])*x[2] * temp - std::sin(x[4])*x[2])*vec[0];
	out[5] = 0 + vec[1] - vec[2];
}

void Crane2DProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
	out[0] = 0;
}

void Crane2DProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
}