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
 * The problem description is based on:
 * Werling, M., Reinisch, P., Gresser, K.: Kombinierte Brems-Ausweich-Assistenz
 * mittels nichtlinearer modellpraediktiver Trajektorienplanung fuer den aktiven
 * Fussgaengerschutz. Tagungsband des 8. Workshop Fahrerassistenzsysteme pp. 68-77 (2012)
 * 
 * Based on GRAMPC:
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * 
 */

#include "vehicle_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

VehicleProblemDescription::VehicleProblemDescription(const std::vector<typeRNum>& pSys, const std::vector<typeRNum>& pCost)
: pSys_(pSys),
  pCost_(pCost)
{
}

void VehicleProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 5;
	*Nu = 2;
	*Np = 0;
	*Nh = 1;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void VehicleProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = std::cos(x[2])*x[4];
	out[1] = std::sin(x[2])*x[4];
	out[2] = (x[3] * x[4]) / (pSys_[0] * (1 + 1 / POW2(pSys_[1])*POW2(x[4])));
	out[3] = u[0];
	out[4] = u[1];
}

void VehicleProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
	out[1] = 0;
	out[2] = (vec[1] * std::cos(x[2]) - vec[0] * std::sin(x[2]))*x[4];
	out[3] = (vec[2] * POW2(pSys_[1])*x[4]) / (pSys_[0] * POW2(pSys_[1]) + pSys_[0] * POW2(x[4]));
	out[4] = vec[0] * std::cos(x[2]) + vec[1] * std::sin(x[2]) + (vec[2] * POW2(pSys_[1])
		/ POW2(POW2(pSys_[1]) + POW2(x[4]))*x[3] * (pSys_[1] - x[4])*(pSys_[1] + x[4])) / pSys_[0];
}

void VehicleProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = vec[3];
	out[1] = vec[4];
}

void VehicleProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}

void VehicleProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    ctypeRNum *udes = param->udes;
    out[0] = pCost_[10] * POW2(u[0] - udes[0])
		+ pCost_[11] * POW2(u[1] - udes[1])
		+ pCost_[0] * POW2(x[0] - xdes[0])
		+ pCost_[1] * POW2(x[1] - xdes[1])
		+ pCost_[2] * POW2(x[2] - xdes[2])
		+ pCost_[3] * POW2(x[3] - xdes[3])
		+ pCost_[4] * POW2(x[4] - xdes[4]);
}

void VehicleProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = 2 * pCost_[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost_[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost_[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost_[3] * (x[3] - xdes[3]);
	out[4] = 2 * pCost_[4] * (x[4] - xdes[4]);
}

void VehicleProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *udes = param->udes;
    out[0] = 2 * pCost_[10] * (u[0] - udes[0]);
	out[1] = 2 * pCost_[11] * (u[1] - udes[1]);
}

void VehicleProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = pCost_[5] * POW2(x[0] - xdes[0])
		+ pCost_[6] * POW2(x[1] - xdes[1])
		+ pCost_[7] * POW2(x[2] - xdes[2])
		+ pCost_[8] * POW2(x[3] - xdes[3])
		+ pCost_[9] * POW2(x[4] - xdes[4]);
}

void VehicleProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    out[0] = 2 * pCost_[5] * (x[0] - xdes[0]);
	out[1] = 2 * pCost_[6] * (x[1] - xdes[1]);
	out[2] = 2 * pCost_[7] * (x[2] - xdes[2]);
	out[3] = 2 * pCost_[8] * (x[3] - xdes[3]);
	out[4] = 2 * pCost_[9] * (x[4] - xdes[4]);
}

void VehicleProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0;
}

void VehicleProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 4 - POW2(-20 + x[0]) - POW2(((typeRNum)-0.5) + x[1]);
}

void VehicleProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = -2 * vec[0] * (-20 + x[0]);
	out[1] = -2 * vec[0] * (((typeRNum)-0.5) + x[1]);
	out[2] = 0;
	out[3] = 0;
	out[4] = 0;
}

void VehicleProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = 0;
	out[1] = 0;
}

void VehicleProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}