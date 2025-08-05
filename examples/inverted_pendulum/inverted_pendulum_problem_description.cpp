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

#include "inverted_pendulum_problem_description.hpp"

using namespace grampc;

InvertedPendulumProblemDescription::InvertedPendulumProblemDescription(const std::vector<typeRNum> &pSys, const std::vector<typeRNum> &pCost, const std::vector<typeRNum> &pCon)
: pSys_(pSys),
  pCost_(pCost),
  pCon_(pCon)
{
}

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void InvertedPendulumProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 4;
    *Nu  = 1;
    *Np  = 0;
    *Nh  = 2;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void InvertedPendulumProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
	out[0] = x[1];
    out[1] = u[0];
    out[2] = x[3];
    out[3] = -(pSys_[0]*x[3]+pSys_[1]*(pSys_[2]*u[0]*std::cos(x[2])+pSys_[2]*pSys_[3]*std::sin(x[2])))/(pSys_[4]+(pSys_[1]*pSys_[1])*pSys_[2]);
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void InvertedPendulumProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
    ctypeRNum  t2 = pSys_[1]*pSys_[1];
    ctypeRNum  t3 = pSys_[2]*t2;
    ctypeRNum  t4 = pSys_[4]+t3;
    ctypeRNum  t5 = 1.0/t4;
    out[0] = 0.0;
    out[1] = vec[0];
    out[2] = -pSys_[1]*pSys_[2]*t5*vec[3]*(pSys_[3]*std::cos(x[2])-u[0]*std::sin(x[2]));
    out[3] = vec[2]-pSys_[0]*t5*vec[3];
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void InvertedPendulumProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
    out[0] = vec[1]-(pSys_[1]*pSys_[2]*vec[3]*std::cos(x[2]))/(pSys_[4]+(pSys_[1]*pSys_[1])*pSys_[2]);
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void InvertedPendulumProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void InvertedPendulumProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = pCost_[8]*POW2(u[0]-udes[0])+pCost_[0]*POW2(x[0]-xdes[0])+pCost_[1]*POW2(x[1]-xdes[1])+pCost_[2]*POW2(x[2]-xdes[2])+pCost_[3]*POW2(x[3]-xdes[3]);
}


/** Gradient dl/dx **/
void InvertedPendulumProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = pCost_[0]*(x[0]-xdes[0])*2.0;
    out[1] = pCost_[1]*(x[1]-xdes[1])*2.0;
    out[2] = pCost_[2]*(x[2]-xdes[2])*2.0;
    out[3] = pCost_[3]*(x[3]-xdes[3])*2.0;
}


/** Gradient dl/du **/
void InvertedPendulumProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    out[0] = pCost_[8]*(u[0]-udes[0])*2.0;
}


/** Gradient dl/dp **/
void InvertedPendulumProblemDescription::dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void InvertedPendulumProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = pCost_[4]*POW2(x[0]-xdes[0])+pCost_[5]*POW2(x[1]-xdes[1])+pCost_[6]*POW2(x[2]-xdes[2])+pCost_[7]*POW2(x[3]-xdes[3]);
}


/** Gradient dV/dx **/
void InvertedPendulumProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = pCost_[4]*(x[0]-xdes[0])*2.0;
    out[1] = pCost_[5]*(x[1]-xdes[1])*2.0;
    out[2] = pCost_[6]*(x[2]-xdes[2])*2.0;
    out[3] = pCost_[7]*(x[3]-xdes[3])*2.0;
}


/** Gradient dV/dp **/
void InvertedPendulumProblemDescription::dVdp(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
}


/** Gradient dV/dT **/
void InvertedPendulumProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void InvertedPendulumProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = x[1]-pCon_[0];
    out[1] = -x[1]-pCon_[0];
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void InvertedPendulumProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
    out[0] = 0.0;
    out[1] = vec[0]-vec[1];
    out[2] = 0.0;
    out[3] = 0.0;
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void InvertedPendulumProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
    out[0] = 0.0;
}

/** Jacobian df/dx in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    ctypeRNum  t2 = pSys_[1]*pSys_[1];
    ctypeRNum  t3 = pSys_[2]*t2;
    ctypeRNum  t4 = pSys_[4]+t3;
    ctypeRNum  t5 = 1.0/t4;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = 1.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = -pSys_[1]*pSys_[2]*t5*(pSys_[3]*std::cos(x[2])-u[0]*std::sin(x[2]));
    out[12] = 0.0;
    out[13] = 0.0;
    out[14] = 1.0;
    out[15] = -pSys_[0]*t5;
}


/** Hessian d(df/dx)/dx in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dfdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = 0.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = 0.0;
    out[12] = 0.0;
    out[13] = 0.0;
    out[14] = 0.0;
    out[15] = 0.0;
    out[16] = 0.0;
    out[17] = 0.0;
    out[18] = 0.0;
    out[19] = 0.0;
    out[20] = 0.0;
    out[21] = 0.0;
    out[22] = 0.0;
    out[23] = 0.0;
    out[24] = 0.0;
    out[25] = 0.0;
    out[26] = 0.0;
    out[27] = 0.0;
    out[28] = 0.0;
    out[29] = 0.0;
    out[30] = 0.0;
    out[31] = 0.0;
    out[32] = 0.0;
    out[33] = 0.0;
    out[34] = 0.0;
    out[35] = 0.0;
    out[36] = 0.0;
    out[37] = 0.0;
    out[38] = 0.0;
    out[39] = 0.0;
    out[40] = 0.0;
    out[41] = 0.0;
    out[42] = 0.0;
    out[43] = (pSys_[1]*pSys_[2]*(u[0]*std::cos(x[2])+pSys_[3]*std::sin(x[2])))/(pSys_[4]+(pSys_[1]*pSys_[1])*pSys_[2]);
    out[44] = 0.0;
    out[45] = 0.0;
    out[46] = 0.0;
    out[47] = 0.0;
    out[48] = 0.0;
    out[49] = 0.0;
    out[50] = 0.0;
    out[51] = 0.0;
    out[52] = 0.0;
    out[53] = 0.0;
    out[54] = 0.0;
    out[55] = 0.0;
    out[56] = 0.0;
    out[57] = 0.0;
    out[58] = 0.0;
    out[59] = 0.0;
    out[60] = 0.0;
    out[61] = 0.0;
    out[62] = 0.0;
    out[63] = 0.0;
}


/** Jacobian d(df/dx)/du in vector form (column-wise **/
void InvertedPendulumProblemDescription::dfdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = 0.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = (pSys_[1]*pSys_[2]*std::sin(x[2]))/(pSys_[4]+(pSys_[1]*pSys_[1])*pSys_[2]);
    out[12] = 0.0;
    out[13] = 0.0;
    out[14] = 0.0;
    out[15] = 0.0;
}


/** Jacobian d(df/dx)/dp in vector form (column-wise **/
void InvertedPendulumProblemDescription::dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian df/dp in vector form (column-wise **/
void InvertedPendulumProblemDescription::dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian d(df/dp)/du in vector form (column-wise **/
void InvertedPendulumProblemDescription::dfdpdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian dh/dx in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dhdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 1.0;
    out[3] = -1.0;
    out[4] = 0.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
}


/** Jacobian dh/du in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dhdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
}

/** Hessian d(dh/dx)/dx in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dhdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = 0.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = 0.0;
    out[12] = 0.0;
    out[13] = 0.0;
    out[14] = 0.0;
    out[15] = 0.0;
    out[16] = 0.0;
    out[17] = 0.0;
    out[18] = 0.0;
    out[19] = 0.0;
    out[20] = 0.0;
    out[21] = 0.0;
    out[22] = 0.0;
    out[23] = 0.0;
    out[24] = 0.0;
    out[25] = 0.0;
    out[26] = 0.0;
    out[27] = 0.0;
    out[28] = 0.0;
    out[29] = 0.0;
    out[30] = 0.0;
    out[31] = 0.0;
}


/** Jacobian d(dh/dx)/du in vector form (column-wise) **/
void InvertedPendulumProblemDescription::dhdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = 0.0;
    out[5] = 0.0;
    out[6] = 0.0;
    out[7] = 0.0;
}