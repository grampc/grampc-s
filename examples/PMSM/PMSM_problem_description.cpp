/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 * This probfct-file describes a modified motor (PMSM) problem from
 * Englert, T., Graichen, K.: Model Predictive Torque Control of PMSMs for
 * High Performance Applications. Control Engineering Practice, Volume 81, 2018, p. 43-54.
 * 
 * Based on GRAMPC:
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * 
 */

#include "PMSM_problem_description.hpp"

#define POW2(a) ((a)*(a))

PMSMProblemDescription::PMSMProblemDescription(const std::vector<typeRNum>& pSys, const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
: pSys_(pSys),
  pCost_(pCost),
  pCon_(pCon)
{
}

void PMSMProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 4;
	*Nu = 2;
	*Np = 3;
	*Nh = 2;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void PMSMProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = (u[0] - pSys_[0] * x[0] + p[1] * x[1] * x[2]) / p[0];
	out[1] = (u[1] - pSys_[0] * x[1] - (p[2] + p[0] * x[0])*x[2]) / p[1];
	out[2] = ((typeRNum)0.5*(-2 * pSys_[4] * pSys_[7] + 3 * POW2(pSys_[4])*
		(p[2] + (p[0] - p[1])*x[0])*x[1] - 2 * pSys_[6] * x[2])) / pSys_[5];
	out[3] = x[2];
}

void PMSMProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = -((vec[0] * pSys_[0]) / p[0]) + ((typeRNum)1.5*vec[2] * POW2(pSys_[4])*(p[0] - p[1])*x[1]) /
		pSys_[5] - (vec[1] * p[0] * x[2]) / p[1];
	out[1] = -((vec[1] * pSys_[0]) / p[1]) + ((typeRNum)1.5*vec[2] * POW2(pSys_[4])*(p[2] + (p[0] - p[1])*x[0])) /
		pSys_[5] + (vec[0] * p[1] * x[2]) / p[0];
	out[2] = vec[3] - (vec[2] * pSys_[6]) / pSys_[5] - (vec[1] * (p[2] + p[0] * x[0])) /
		p[1] + (vec[0] * p[1] * x[1]) / p[0];
	out[3] = 0;
}

void PMSMProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
	out[0] = vec[0] / p[0];
	out[1] = vec[1] / p[1];
}

void PMSMProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
	out[0] = - vec[0] * (u[0] - pSys_[0] * x[0] + p[1] * x[1] * x[2]) / p[0] / p[0]
			 - vec[1] * x[0]*x[2] / p[1]
			 + vec[2] * 1.5 * POW2(pSys_[4]) * x[0]*x[1] / pSys_[5];
	out[1] =   vec[0] * x[1] * x[2] / p[0]
			 - vec[1] * (u[1] - pSys_[0] * x[1] - (p[2] + p[0] * x[0])*x[2]) / p[1] / p[1]
			 - vec[2] * 1.5 * POW2(pSys_[4]) * x[0]*x[1] / pSys_[5];
	out[2] = - vec[1] * x[2] / p[1]
			 + vec[2] * 1.5 * POW2(pSys_[4]) * x[1] / pSys_[5];
}


void PMSMProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    ctypeRNum *udes = param->udes;
    out[0] = pCost_[0] * POW2(x[0] - xdes[0])
		+ pCost_[1] * POW2(x[1] - xdes[1])
		+ pCost_[2] * POW2(x[2] - xdes[2])
		+ pCost_[3] * POW2(x[3] - xdes[3])
		+ pCost_[4] * POW2(u[0] - udes[0])
		+ pCost_[5] * POW2(u[1] - udes[1]);
}

void PMSMProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = 2 * pCost_[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost_[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost_[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost_[3] * (x[3] - xdes[3]);
}

void PMSMProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *udes = param->udes;
    out[0] = 2 * pCost_[4] * (u[0] - udes[0]);
	out[1] = 2 * pCost_[5] * (u[1] - udes[1]);
}

void PMSMProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}

void PMSMProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}

void PMSMProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}

void PMSMProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = ((POW2(u[0]) + POW2(u[1])) - pCon_[0]) / pCon_[0];
	out[1] = ((POW2(x[0]) + POW2(x[1])) - pCon_[1]) / pCon_[1];
}

void PMSMProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = (2 * vec[1] * x[0]) / pCon_[1];
	out[1] = (2 * vec[1] * x[1]) / pCon_[1];
	out[2] = 0;
	out[3] = 0;
}

void PMSMProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = (2 * vec[0] * u[0]) / pCon_[0];
	out[1] = (2 * vec[0] * u[1]) / pCon_[0];
}

void PMSMProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
	out[1] = 0.0;
	out[2] = 0.0;
}