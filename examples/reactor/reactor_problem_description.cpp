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
 * This probfct-file describes a simplified nonlinear chemical reactor from
 * R. Rothfuss, J. Rudolph, M. Zeitz: Flatness based control of a nonlinear chemical reactor model. 
 * Automatica, Volume 32, Issue 10, 1996, p. 1433-1439.
 * 
 */

#include "reactor_problem_description.hpp"

ReactorProblemDescription::ReactorProblemDescription(const std::vector<typeRNum>& pSys,
                                                     const std::vector<typeRNum>& pCost,
                                                     const std::vector<typeRNum>& pCon)
    : pSys_(pSys), pCost_(pCost), pCon_(pCon)
{
}

void ReactorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 2;
    *Nu = 1;
    *Np = 0;
    *Ng = 0;
    *Nh = 1;
    *NgT = 0;
    *NhT = 0;
}

void ReactorProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = -pSys_[0] * x[0] - pSys_[2] * x[0] * x[0] + (1 - x[0]) * u[0];
    out[1] = pSys_[0] * x[0] - pSys_[1] * x[1] - x[1] * u[0];
}

void ReactorProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
{
    out[0] = (-pSys_[0] - pSys_[2] * 2 * x[0] - u[0]) * adj[0] + pSys_[0] * adj[1];
    out[1] = (-pSys_[1] - u[0]) * adj[1];
}

void ReactorProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p)
{
    out[0] = (1 - x[0]) * adj[0] + (-x[1]) * adj[1];
}

void ReactorProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = pCost_[2] * (x[0] - xdes[0]) * (x[0] - xdes[0]) +
             pCost_[3] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +
             pCost_[4] * (u[0] - udes[0]) * (u[0] - udes[0]);
}

void ReactorProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = 2 * pCost_[2] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[3] * (x[1] - xdes[1]);
}

void ReactorProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = 2 * pCost_[4] * (u[0] - udes[0]);
}

void ReactorProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = pCost_[0] * (x[0] - xdes[0]) * (x[0] - xdes[0]) +
             pCost_[1] * (x[1] - xdes[1]) * (x[1] - xdes[1]);
}

void ReactorProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 2 * pCost_[0] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[1] * (x[1] - xdes[1]);
}

void ReactorProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = x[1] - pCon_[0];
}

void ReactorProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
    out[0] = 0.0;
    out[1] = vec[0];
}

void ReactorProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
    out[0] = 0.0;
}
