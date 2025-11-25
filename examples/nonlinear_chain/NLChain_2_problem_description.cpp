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
* This probfct-file describes the nonlinear chain problem with 2 chain elements from
* Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
* Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
* Process Control 22(3), 540-550 (2012)
*
*/

#include "NLChain_2_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

using namespace grampc;


Chain_2_ProblemDescription::Chain_2_ProblemDescription(const std::vector<typeRNum>& pCost)
: pCost_(pCost)
{
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void Chain_2_ProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 9;
    *Nu  = 3;
    *Np  = 0;
    *Nh  = 0;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void Chain_2_ProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum  t2 = x[0]*x[0];
    ctypeRNum  t3 = x[1]*x[1];
    ctypeRNum  t4 = x[2]*x[2];
    ctypeRNum  t5 = -x[6];
    ctypeRNum  t6 = -x[7];
    ctypeRNum  t7 = -x[8];
    ctypeRNum  t8 = t5+x[0];
    ctypeRNum  t9 = t6+x[1];
    ctypeRNum  t10 = t7+x[2];
    ctypeRNum  t14 = t2+t3+t4;
    ctypeRNum  t11 = t8*t8;
    ctypeRNum  t12 = t9*t9;
    ctypeRNum  t13 = t10*t10;
    ctypeRNum  t15 = 1.0/sqrt(t14);
    ctypeRNum  t16 = t15*(1.1E+1/4.0E+2);
    ctypeRNum  t17 = t11+t12+t13;
    ctypeRNum  t18 = t16-1.0/1.0E+1;
    ctypeRNum  t19 = 1.0/sqrt(t17);
    ctypeRNum  t20 = t19*(1.1E+1/4.0E+2);
    ctypeRNum  t21 = t20-1.0/1.0E+1;
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = t8*t21*(1.0E+2/3.0)+t18*x[0]*(1.0E+2/3.0);
    out[4] = t9*t21*(1.0E+2/3.0)+t18*x[1]*(1.0E+2/3.0);
    out[5] = t10*t21*(1.0E+2/3.0)+t18*x[2]*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[6] = u[0];
    out[7] = u[1];
    out[8] = u[2];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Chain_2_ProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    ctypeRNum  t2 = x[0]*x[0];
    ctypeRNum  t3 = x[1]*x[1];
    ctypeRNum  t4 = x[2]*x[2];
    ctypeRNum  t5 = -x[6];
    ctypeRNum  t6 = -x[7];
    ctypeRNum  t7 = -x[8];
    ctypeRNum  t8 = vec[3]*x[0]*1.1E+1;
    ctypeRNum  t9 = vec[4]*x[1]*1.1E+1;
    ctypeRNum  t10 = vec[5]*x[2]*1.1E+1;
    ctypeRNum  t11 = vec[3]*x[6]*1.1E+1;
    ctypeRNum  t12 = vec[4]*x[7]*1.1E+1;
    ctypeRNum  t13 = vec[5]*x[8]*1.1E+1;
    ctypeRNum  t14 = -t11;
    ctypeRNum  t15 = -t12;
    ctypeRNum  t16 = -t13;
    ctypeRNum  t17 = t5+x[0];
    ctypeRNum  t18 = t6+x[1];
    ctypeRNum  t19 = t7+x[2];
    ctypeRNum  t23 = t2+t3+t4;
    ctypeRNum  t20 = t17*t17;
    ctypeRNum  t21 = t18*t18;
    ctypeRNum  t22 = t19*t19;
    ctypeRNum  t24 = 1.0/sqrt(t23);
    ctypeRNum  t25 = t24*t24*t24;
    ctypeRNum  t26 = t24*(1.1E+1/1.2E+1);
    ctypeRNum  t31 = t20+t21+t22;
    ctypeRNum  t27 = -t26;
    ctypeRNum  t28 = t25*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t29 = t25*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t30 = t25*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t32 = pow(t31,3.0/2.0);
    ctypeRNum  t33 = 1.0/sqrt(t31);
    ctypeRNum  t34 = 1.0/t32;
    ctypeRNum  t35 = t33*(1.1E+1/1.2E+1);
    ctypeRNum  t36 = -t35;
    ctypeRNum  t37 = t17*t18*t34*(1.1E+1/1.2E+1);
    ctypeRNum  t38 = t17*t19*t34*(1.1E+1/1.2E+1);
    ctypeRNum  t39 = t18*t19*t34*(1.1E+1/1.2E+1);
    ctypeRNum  t40 = t28+t37;
    ctypeRNum  t41 = t29+t38;
    ctypeRNum  t42 = t30+t39;
    out[0] = -t40*vec[4]-t41*vec[5]-vec[3]*(t27+t36+t2*t25*(1.1E+1/1.2E+1)+t17*t34*(x[0]*2.0-x[6]*2.0)*(1.1E+1/2.4E+1)+2.0E+1/3.0);
    out[1] = -t40*vec[3]-t42*vec[5]-vec[4]*(t27+t36+t3*t25*(1.1E+1/1.2E+1)+t18*t34*(x[1]*2.0-x[7]*2.0)*(1.1E+1/2.4E+1)+2.0E+1/3.0);
    out[2] = -t41*vec[3]-t42*vec[4]-vec[5]*(t27+t36+t4*t25*(1.1E+1/1.2E+1)+t19*t34*(x[2]*2.0-x[8]*2.0)*(1.1E+1/2.4E+1)+2.0E+1/3.0);
    out[3] = vec[0];
    out[4] = vec[1];
    out[5] = vec[2];
    out[6] = (t34*(t2*vec[3]*1.1E+1+t32*vec[3]*4.0E+1+t9*x[0]+t10*x[0]-vec[3]*(t2+t21+t22)*1.1E+1-vec[4]*x[0]*x[7]*1.1E+1-vec[5]*x[0]*x[8]*1.1E+1))/1.2E+1-(t34*x[6]*(t9+t10+t15+t16))/1.2E+1;
    out[7] = (t34*(t3*vec[4]*1.1E+1+t32*vec[4]*4.0E+1+t8*x[1]+t10*x[1]-vec[4]*(t3+t20+t22)*1.1E+1-vec[3]*x[1]*x[6]*1.1E+1-vec[5]*x[1]*x[8]*1.1E+1))/1.2E+1-(t34*x[7]*(t8+t10+t14+t16))/1.2E+1;
    out[8] = (t34*(t4*vec[5]*1.1E+1+t32*vec[5]*4.0E+1+t8*x[2]+t9*x[2]-vec[5]*(t4+t20+t21)*1.1E+1-vec[3]*x[2]*x[6]*1.1E+1-vec[4]*x[2]*x[7]*1.1E+1))/1.2E+1-(t34*x[8]*(t8+t9+t14+t15))/1.2E+1;
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Chain_2_ProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = vec[6];
    out[1] = vec[7];
    out[2] = vec[8];
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void Chain_2_ProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void Chain_2_ProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*(u[0]*u[0])+pCost_[2]*(u[1]*u[1])+pCost_[2]*(u[2]*u[2])+pCost_[1]*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+pCost_[0]*(POW2(x[6]-1.0)+x[7]*x[7]+x[8]*x[8]);
}


/** Gradient dl/dx **/
void Chain_2_ProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = pCost_[1]*x[3]*2.0;
    out[4] = pCost_[1]*x[4]*2.0;
    out[5] = pCost_[1]*x[5]*2.0;
    out[6] = pCost_[0]*(x[6]-1.0)*2.0;
    out[7] = pCost_[0]*x[7]*2.0;
    out[8] = pCost_[0]*x[8]*2.0;
}


/** Gradient dl/du **/
void Chain_2_ProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*u[0]*2.0;
    out[1] = pCost_[2]*u[1]*2.0;
    out[2] = pCost_[2]*u[2]*2.0;
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Chain_2_ProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Gradient dV/dx **/
void Chain_2_ProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
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
}


/** Gradient dV/dT **/
void Chain_2_ProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void Chain_2_ProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Chain_2_ProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Chain_2_ProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void Chain_2_ProblemDescription::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void Chain_2_ProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void Chain_2_ProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian df/dx in vector form (column-wise) **/
void Chain_2_ProblemDescription::dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum  t2 = x[0]*2.0;
    ctypeRNum  t3 = x[1]*2.0;
    ctypeRNum  t4 = x[2]*2.0;
    ctypeRNum  t5 = x[6]*2.0;
    ctypeRNum  t6 = x[7]*2.0;
    ctypeRNum  t7 = x[8]*2.0;
    ctypeRNum  t8 = x[0]*x[0];
    ctypeRNum  t9 = x[1]*x[1];
    ctypeRNum  t10 = x[2]*x[2];
    ctypeRNum  t11 = -x[6];
    ctypeRNum  t13 = -x[7];
    ctypeRNum  t15 = -x[8];
    ctypeRNum  t12 = -t5;
    ctypeRNum  t14 = -t6;
    ctypeRNum  t16 = -t7;
    ctypeRNum  t17 = t11+x[0];
    ctypeRNum  t18 = t13+x[1];
    ctypeRNum  t19 = t15+x[2];
    ctypeRNum  t26 = t8+t9+t10;
    ctypeRNum  t20 = t2+t12;
    ctypeRNum  t21 = t3+t14;
    ctypeRNum  t22 = t4+t16;
    ctypeRNum  t23 = t17*t17;
    ctypeRNum  t24 = t18*t18;
    ctypeRNum  t25 = t19*t19;
    ctypeRNum  t27 = 1.0/sqrt(t26);
    ctypeRNum  t28 = t27*t27*t27;
    ctypeRNum  t29 = t27*(1.1E+1/1.2E+1);
    ctypeRNum  t36 = t23+t24+t25;
    ctypeRNum  t30 = t28*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t31 = t28*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t32 = t28*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t37 = 1.0/sqrt(t36);
    ctypeRNum  t33 = -t30;
    ctypeRNum  t34 = -t31;
    ctypeRNum  t35 = -t32;
    ctypeRNum  t38 = t37*t37*t37;
    ctypeRNum  t39 = t37*(1.1E+1/1.2E+1);
    ctypeRNum  t40 = -t39;
    ctypeRNum  t41 = t17*t18*t38*(1.1E+1/1.2E+1);
    ctypeRNum  t42 = t17*t19*t38*(1.1E+1/1.2E+1);
    ctypeRNum  t43 = t18*t19*t38*(1.1E+1/1.2E+1);
    ctypeRNum  t47 = t17*t20*t38*(1.1E+1/2.4E+1);
    ctypeRNum  t48 = t18*t21*t38*(1.1E+1/2.4E+1);
    ctypeRNum  t49 = t19*t22*t38*(1.1E+1/2.4E+1);
    ctypeRNum  t44 = -t41;
    ctypeRNum  t45 = -t42;
    ctypeRNum  t46 = -t43;
    ctypeRNum  t50 = t33+t44;
    ctypeRNum  t51 = t34+t45;
    ctypeRNum  t52 = t35+t46;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = t29+t39-t47-t8*t28*(1.1E+1/1.2E+1)-2.0E+1/3.0;
    out[4] = t50;
    out[5] = t51;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = 0.0;
    out[12] = t50;
    out[13] = t29+t39-t48-t9*t28*(1.1E+1/1.2E+1)-2.0E+1/3.0;
    out[14] = t52;
    out[15] = 0.0;
    out[16] = 0.0;
    out[17] = 0.0;
    out[18] = 0.0;
    out[19] = 0.0;
    out[20] = 0.0;
    out[21] = t51;
    out[22] = t52;
    out[23] = t29+t39-t49-t10*t28*(1.1E+1/1.2E+1)-2.0E+1/3.0;
    out[24] = 0.0;
    out[25] = 0.0;
    out[26] = 0.0;
    out[27] = 1.0;
    out[28] = 0.0;
    out[29] = 0.0;
    out[30] = 0.0;
    out[31] = 0.0;
    out[32] = 0.0;
    out[33] = 0.0;
    out[34] = 0.0;
    out[35] = 0.0;
    out[36] = 0.0;
    out[37] = 1.0;
    out[38] = 0.0;
    out[39] = 0.0;
    out[40] = 0.0;
    out[41] = 0.0;
    out[42] = 0.0;
    out[43] = 0.0;
    out[44] = 0.0;
    out[45] = 0.0;
    out[46] = 0.0;
    out[47] = 1.0;
    out[48] = 0.0;
    out[49] = 0.0;
    out[50] = 0.0;
    out[51] = 0.0;
    out[52] = 0.0;
    out[53] = 0.0;
    out[54] = 0.0;
    out[55] = 0.0;
    out[56] = 0.0;
    out[57] = t40+t47+1.0E+1/3.0;
    out[58] = t41;
    out[59] = t42;
    out[60] = 0.0;
    out[61] = 0.0;
    out[62] = 0.0;
    out[63] = 0.0;
    out[64] = 0.0;
    out[65] = 0.0;
    out[66] = t41;
    out[67] = t40+t48+1.0E+1/3.0;
    out[68] = t43;
    out[69] = 0.0;
    out[70] = 0.0;
    out[71] = 0.0;
    out[72] = 0.0;
    out[73] = 0.0;
    out[74] = 0.0;
    out[75] = t42;
    out[76] = t43;
    out[77] = t40+t49+1.0E+1/3.0;
    out[78] = 0.0;
    out[79] = 0.0;
    out[80] = 0.0;
}


/** Hessian d(df/dx)/dx in vector form (column-wise) **/
void Chain_2_ProblemDescription::dfdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    ctypeRNum  t2 = x[0]*2.0;
    ctypeRNum  t3 = x[1]*2.0;
    ctypeRNum  t4 = x[2]*2.0;
    ctypeRNum  t5 = x[6]*2.0;
    ctypeRNum  t6 = x[7]*2.0;
    ctypeRNum  t7 = x[8]*2.0;
    ctypeRNum  t8 = x[0]*x[0];
    ctypeRNum  t9 = x[1]*x[1];
    ctypeRNum  t10 = x[2]*x[2];
    ctypeRNum  t11 = x[6]*x[6];
    ctypeRNum  t12 = x[7]*x[7];
    ctypeRNum  t13 = x[8]*x[8];
    ctypeRNum  t15 = x[0]*x[6]*4.0;
    ctypeRNum  t17 = x[1]*x[7]*4.0;
    ctypeRNum  t19 = x[2]*x[8]*4.0;
    ctypeRNum  t20 = -x[6];
    ctypeRNum  t22 = -x[7];
    ctypeRNum  t24 = -x[8];
    ctypeRNum  t32 = x[0]*x[6]*-2.0;
    ctypeRNum  t33 = x[1]*x[7]*-2.0;
    ctypeRNum  t34 = x[2]*x[8]*-2.0;
    ctypeRNum  t21 = -t5;
    ctypeRNum  t23 = -t6;
    ctypeRNum  t25 = -t7;
    ctypeRNum  t35 = t8*-2.0;
    ctypeRNum  t36 = t9*-2.0;
    ctypeRNum  t37 = t10*-2.0;
    ctypeRNum  t38 = t11*-2.0;
    ctypeRNum  t39 = t12*-2.0;
    ctypeRNum  t40 = t13*-2.0;
    ctypeRNum  t41 = t20+x[0];
    ctypeRNum  t42 = t22+x[1];
    ctypeRNum  t43 = t24+x[2];
    ctypeRNum  t50 = t8+t9+t10;
    ctypeRNum  t60 = t8+t9+t11+t12+t32+t33;
    ctypeRNum  t61 = t8+t10+t11+t13+t32+t34;
    ctypeRNum  t62 = t9+t10+t12+t13+t33+t34;
    ctypeRNum  t44 = t2+t21;
    ctypeRNum  t45 = t3+t23;
    ctypeRNum  t46 = t4+t25;
    ctypeRNum  t47 = t41*t41;
    ctypeRNum  t48 = t42*t42;
    ctypeRNum  t49 = t43*t43;
    ctypeRNum  t51 = 1.0/pow(t50,3.0/2.0);
    ctypeRNum  t52 = 1.0/pow(t50,5.0/2.0);
    ctypeRNum  t72 = t19+t37+t40+t60;
    ctypeRNum  t73 = t17+t36+t39+t61;
    ctypeRNum  t74 = t15+t35+t38+t62;
    ctypeRNum  t53 = t51*x[0]*(1.1E+1/1.2E+1);
    ctypeRNum  t54 = t51*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t55 = t51*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t56 = t52*x[0]*x[1]*x[2]*(1.1E+1/4.0);
    ctypeRNum  t63 = t9*t52*x[0]*(1.1E+1/4.0);
    ctypeRNum  t64 = t8*t52*x[1]*(1.1E+1/4.0);
    ctypeRNum  t65 = t10*t52*x[0]*(1.1E+1/4.0);
    ctypeRNum  t66 = t8*t52*x[2]*(1.1E+1/4.0);
    ctypeRNum  t67 = t10*t52*x[1]*(1.1E+1/4.0);
    ctypeRNum  t68 = t9*t52*x[2]*(1.1E+1/4.0);
    ctypeRNum  t69 = t47+t48+t49;
    ctypeRNum  t57 = -t53;
    ctypeRNum  t58 = -t54;
    ctypeRNum  t59 = -t55;
    ctypeRNum  t70 = 1.0/pow(t69,3.0/2.0);
    ctypeRNum  t71 = 1.0/pow(t69,5.0/2.0);
    ctypeRNum  t75 = t41*t70*(1.1E+1/1.2E+1);
    ctypeRNum  t76 = t42*t70*(1.1E+1/1.2E+1);
    ctypeRNum  t77 = t43*t70*(1.1E+1/1.2E+1);
    ctypeRNum  t81 = t44*t70*(1.1E+1/2.4E+1);
    ctypeRNum  t82 = t45*t70*(1.1E+1/2.4E+1);
    ctypeRNum  t83 = t46*t70*(1.1E+1/2.4E+1);
    ctypeRNum  t87 = t41*t42*t43*t71*(1.1E+1/4.0);
    ctypeRNum  t89 = t41*t42*t44*t71*(1.1E+1/8.0);
    ctypeRNum  t90 = t41*t43*t44*t71*(1.1E+1/8.0);
    ctypeRNum  t91 = t41*t42*t45*t71*(1.1E+1/8.0);
    ctypeRNum  t92 = t42*t43*t45*t71*(1.1E+1/8.0);
    ctypeRNum  t93 = t41*t43*t46*t71*(1.1E+1/8.0);
    ctypeRNum  t94 = t42*t43*t46*t71*(1.1E+1/8.0);
    ctypeRNum  t95 = t43*t60*t71*(1.1E+1/4.0);
    ctypeRNum  t96 = t42*t61*t71*(1.1E+1/4.0);
    ctypeRNum  t97 = t41*t62*t71*(1.1E+1/4.0);
    ctypeRNum  t99 = t41*t71*t72*(1.1E+1/1.2E+1);
    ctypeRNum  t100 = t41*t71*t73*(1.1E+1/1.2E+1);
    ctypeRNum  t101 = t42*t71*t72*(1.1E+1/1.2E+1);
    ctypeRNum  t102 = t42*t71*t74*(1.1E+1/1.2E+1);
    ctypeRNum  t103 = t43*t71*t73*(1.1E+1/1.2E+1);
    ctypeRNum  t104 = t43*t71*t74*(1.1E+1/1.2E+1);
    ctypeRNum  t78 = -t75;
    ctypeRNum  t79 = -t76;
    ctypeRNum  t80 = -t77;
    ctypeRNum  t84 = -t81;
    ctypeRNum  t85 = -t82;
    ctypeRNum  t86 = -t83;
    ctypeRNum  t88 = -t87;
    ctypeRNum  t98 = t56+t87;
    ctypeRNum  t105 = -t99;
    ctypeRNum  t106 = -t100;
    ctypeRNum  t107 = -t101;
    ctypeRNum  t108 = -t102;
    ctypeRNum  t109 = -t103;
    ctypeRNum  t110 = -t104;
    ctypeRNum  t111 = t57+t63+t78+t91;
    ctypeRNum  t112 = t58+t64+t79+t89;
    ctypeRNum  t113 = t57+t65+t78+t93;
    ctypeRNum  t114 = t59+t66+t80+t90;
    ctypeRNum  t115 = t58+t67+t79+t94;
    ctypeRNum  t116 = t59+t68+t80+t92;
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = t78-t44*t70*(1.1E+1/1.2E+1)-t51*x[0]*(1.1E+1/4.0)+(t41*t41*t41)*t71*(1.1E+1/4.0)+t52*(x[0]*x[0]*x[0])*(1.1E+1/4.0);
    out[4] = t112;
    out[5] = t114;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = 0.0;
    out[10] = 0.0;
    out[11] = 0.0;
    out[12] = t112;
    out[13] = t57+t63+t84+t42*t44*t45*t71*(1.1E+1/1.6E+1);
    out[14] = t98;
    out[15] = 0.0;
    out[16] = 0.0;
    out[17] = 0.0;
    out[18] = 0.0;
    out[19] = 0.0;
    out[20] = 0.0;
    out[21] = t114;
    out[22] = t98;
    out[23] = t57+t65+t84+t43*t44*t46*t71*(1.1E+1/1.6E+1);
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
    out[43] = 0.0;
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
    out[57] = t97;
    out[58] = t102;
    out[59] = t104;
    out[60] = 0.0;
    out[61] = 0.0;
    out[62] = 0.0;
    out[63] = 0.0;
    out[64] = 0.0;
    out[65] = 0.0;
    out[66] = t102;
    out[67] = t100;
    out[68] = t88;
    out[69] = 0.0;
    out[70] = 0.0;
    out[71] = 0.0;
    out[72] = 0.0;
    out[73] = 0.0;
    out[74] = 0.0;
    out[75] = t104;
    out[76] = t88;
    out[77] = t99;
    out[78] = 0.0;
    out[79] = 0.0;
    out[80] = 0.0;
    out[81] = 0.0;
    out[82] = 0.0;
    out[83] = 0.0;
    out[84] = t58+t64+t85+t41*t44*t45*t71*(1.1E+1/1.6E+1);
    out[85] = t111;
    out[86] = t98;
    out[87] = 0.0;
    out[88] = 0.0;
    out[89] = 0.0;
    out[90] = 0.0;
    out[91] = 0.0;
    out[92] = 0.0;
    out[93] = t111;
    out[94] = t79-t45*t70*(1.1E+1/1.2E+1)-t51*x[1]*(1.1E+1/4.0)+(t42*t42*t42)*t71*(1.1E+1/4.0)+t52*(x[1]*x[1]*x[1])*(1.1E+1/4.0);
    out[95] = t116;
    out[96] = 0.0;
    out[97] = 0.0;
    out[98] = 0.0;
    out[99] = 0.0;
    out[100] = 0.0;
    out[101] = 0.0;
    out[102] = t98;
    out[103] = t116;
    out[104] = t58+t67+t85+t43*t45*t46*t71*(1.1E+1/1.6E+1);
    out[105] = 0.0;
    out[106] = 0.0;
    out[107] = 0.0;
    out[108] = 0.0;
    out[109] = 0.0;
    out[110] = 0.0;
    out[111] = 0.0;
    out[112] = 0.0;
    out[113] = 0.0;
    out[114] = 0.0;
    out[115] = 0.0;
    out[116] = 0.0;
    out[117] = 0.0;
    out[118] = 0.0;
    out[119] = 0.0;
    out[120] = 0.0;
    out[121] = 0.0;
    out[122] = 0.0;
    out[123] = 0.0;
    out[124] = 0.0;
    out[125] = 0.0;
    out[126] = 0.0;
    out[127] = 0.0;
    out[128] = 0.0;
    out[129] = 0.0;
    out[130] = 0.0;
    out[131] = 0.0;
    out[132] = 0.0;
    out[133] = 0.0;
    out[134] = 0.0;
    out[135] = 0.0;
    out[136] = 0.0;
    out[137] = 0.0;
    out[138] = t102;
    out[139] = t100;
    out[140] = t88;
    out[141] = 0.0;
    out[142] = 0.0;
    out[143] = 0.0;
    out[144] = 0.0;
    out[145] = 0.0;
    out[146] = 0.0;
    out[147] = t100;
    out[148] = t96;
    out[149] = t103;
    out[150] = 0.0;
    out[151] = 0.0;
    out[152] = 0.0;
    out[153] = 0.0;
    out[154] = 0.0;
    out[155] = 0.0;
    out[156] = t88;
    out[157] = t103;
    out[158] = t101;
    out[159] = 0.0;
    out[160] = 0.0;
    out[161] = 0.0;
    out[162] = 0.0;
    out[163] = 0.0;
    out[164] = 0.0;
    out[165] = t59+t66+t86+t41*t44*t46*t71*(1.1E+1/1.6E+1);
    out[166] = t98;
    out[167] = t113;
    out[168] = 0.0;
    out[169] = 0.0;
    out[170] = 0.0;
    out[171] = 0.0;
    out[172] = 0.0;
    out[173] = 0.0;
    out[174] = t98;
    out[175] = t59+t68+t86+t42*t45*t46*t71*(1.1E+1/1.6E+1);
    out[176] = t115;
    out[177] = 0.0;
    out[178] = 0.0;
    out[179] = 0.0;
    out[180] = 0.0;
    out[181] = 0.0;
    out[182] = 0.0;
    out[183] = t113;
    out[184] = t115;
    out[185] = t80-t46*t70*(1.1E+1/1.2E+1)-t51*x[2]*(1.1E+1/4.0)+(t43*t43*t43)*t71*(1.1E+1/4.0)+t52*(x[2]*x[2]*x[2])*(1.1E+1/4.0);
    out[186] = 0.0;
    out[187] = 0.0;
    out[188] = 0.0;
    out[189] = 0.0;
    out[190] = 0.0;
    out[191] = 0.0;
    out[192] = 0.0;
    out[193] = 0.0;
    out[194] = 0.0;
    out[195] = 0.0;
    out[196] = 0.0;
    out[197] = 0.0;
    out[198] = 0.0;
    out[199] = 0.0;
    out[200] = 0.0;
    out[201] = 0.0;
    out[202] = 0.0;
    out[203] = 0.0;
    out[204] = 0.0;
    out[205] = 0.0;
    out[206] = 0.0;
    out[207] = 0.0;
    out[208] = 0.0;
    out[209] = 0.0;
    out[210] = 0.0;
    out[211] = 0.0;
    out[212] = 0.0;
    out[213] = 0.0;
    out[214] = 0.0;
    out[215] = 0.0;
    out[216] = 0.0;
    out[217] = 0.0;
    out[218] = 0.0;
    out[219] = t104;
    out[220] = t88;
    out[221] = t99;
    out[222] = 0.0;
    out[223] = 0.0;
    out[224] = 0.0;
    out[225] = 0.0;
    out[226] = 0.0;
    out[227] = 0.0;
    out[228] = t88;
    out[229] = t103;
    out[230] = t101;
    out[231] = 0.0;
    out[232] = 0.0;
    out[233] = 0.0;
    out[234] = 0.0;
    out[235] = 0.0;
    out[236] = 0.0;
    out[237] = t99;
    out[238] = t101;
    out[239] = t95;
    out[240] = 0.0;
    out[241] = 0.0;
    out[242] = 0.0;
    out[243] = 0.0;
    out[244] = 0.0;
    out[245] = 0.0;
    out[246] = 0.0;
    out[247] = 0.0;
    out[248] = 0.0;
    out[249] = 0.0;
    out[250] = 0.0;
    out[251] = 0.0;
    out[252] = 0.0;
    out[253] = 0.0;
    out[254] = 0.0;
    out[255] = 0.0;
    out[256] = 0.0;
    out[257] = 0.0;
    out[258] = 0.0;
    out[259] = 0.0;
    out[260] = 0.0;
    out[261] = 0.0;
    out[262] = 0.0;
    out[263] = 0.0;
    out[264] = 0.0;
    out[265] = 0.0;
    out[266] = 0.0;
    out[267] = 0.0;
    out[268] = 0.0;
    out[269] = 0.0;
    out[270] = 0.0;
    out[271] = 0.0;
    out[272] = 0.0;
    out[273] = 0.0;
    out[274] = 0.0;
    out[275] = 0.0;
    out[276] = 0.0;
    out[277] = 0.0;
    out[278] = 0.0;
    out[279] = 0.0;
    out[280] = 0.0;
    out[281] = 0.0;
    out[282] = 0.0;
    out[283] = 0.0;
    out[284] = 0.0;
    out[285] = 0.0;
    out[286] = 0.0;
    out[287] = 0.0;
    out[288] = 0.0;
    out[289] = 0.0;
    out[290] = 0.0;
    out[291] = 0.0;
    out[292] = 0.0;
    out[293] = 0.0;
    out[294] = 0.0;
    out[295] = 0.0;
    out[296] = 0.0;
    out[297] = 0.0;
    out[298] = 0.0;
    out[299] = 0.0;
    out[300] = 0.0;
    out[301] = 0.0;
    out[302] = 0.0;
    out[303] = 0.0;
    out[304] = 0.0;
    out[305] = 0.0;
    out[306] = 0.0;
    out[307] = 0.0;
    out[308] = 0.0;
    out[309] = 0.0;
    out[310] = 0.0;
    out[311] = 0.0;
    out[312] = 0.0;
    out[313] = 0.0;
    out[314] = 0.0;
    out[315] = 0.0;
    out[316] = 0.0;
    out[317] = 0.0;
    out[318] = 0.0;
    out[319] = 0.0;
    out[320] = 0.0;
    out[321] = 0.0;
    out[322] = 0.0;
    out[323] = 0.0;
    out[324] = 0.0;
    out[325] = 0.0;
    out[326] = 0.0;
    out[327] = 0.0;
    out[328] = 0.0;
    out[329] = 0.0;
    out[330] = 0.0;
    out[331] = 0.0;
    out[332] = 0.0;
    out[333] = 0.0;
    out[334] = 0.0;
    out[335] = 0.0;
    out[336] = 0.0;
    out[337] = 0.0;
    out[338] = 0.0;
    out[339] = 0.0;
    out[340] = 0.0;
    out[341] = 0.0;
    out[342] = 0.0;
    out[343] = 0.0;
    out[344] = 0.0;
    out[345] = 0.0;
    out[346] = 0.0;
    out[347] = 0.0;
    out[348] = 0.0;
    out[349] = 0.0;
    out[350] = 0.0;
    out[351] = 0.0;
    out[352] = 0.0;
    out[353] = 0.0;
    out[354] = 0.0;
    out[355] = 0.0;
    out[356] = 0.0;
    out[357] = 0.0;
    out[358] = 0.0;
    out[359] = 0.0;
    out[360] = 0.0;
    out[361] = 0.0;
    out[362] = 0.0;
    out[363] = 0.0;
    out[364] = 0.0;
    out[365] = 0.0;
    out[366] = 0.0;
    out[367] = 0.0;
    out[368] = 0.0;
    out[369] = 0.0;
    out[370] = 0.0;
    out[371] = 0.0;
    out[372] = 0.0;
    out[373] = 0.0;
    out[374] = 0.0;
    out[375] = 0.0;
    out[376] = 0.0;
    out[377] = 0.0;
    out[378] = 0.0;
    out[379] = 0.0;
    out[380] = 0.0;
    out[381] = 0.0;
    out[382] = 0.0;
    out[383] = 0.0;
    out[384] = 0.0;
    out[385] = 0.0;
    out[386] = 0.0;
    out[387] = 0.0;
    out[388] = 0.0;
    out[389] = 0.0;
    out[390] = 0.0;
    out[391] = 0.0;
    out[392] = 0.0;
    out[393] = 0.0;
    out[394] = 0.0;
    out[395] = 0.0;
    out[396] = 0.0;
    out[397] = 0.0;
    out[398] = 0.0;
    out[399] = 0.0;
    out[400] = 0.0;
    out[401] = 0.0;
    out[402] = 0.0;
    out[403] = 0.0;
    out[404] = 0.0;
    out[405] = 0.0;
    out[406] = 0.0;
    out[407] = 0.0;
    out[408] = 0.0;
    out[409] = 0.0;
    out[410] = 0.0;
    out[411] = 0.0;
    out[412] = 0.0;
    out[413] = 0.0;
    out[414] = 0.0;
    out[415] = 0.0;
    out[416] = 0.0;
    out[417] = 0.0;
    out[418] = 0.0;
    out[419] = 0.0;
    out[420] = 0.0;
    out[421] = 0.0;
    out[422] = 0.0;
    out[423] = 0.0;
    out[424] = 0.0;
    out[425] = 0.0;
    out[426] = 0.0;
    out[427] = 0.0;
    out[428] = 0.0;
    out[429] = 0.0;
    out[430] = 0.0;
    out[431] = 0.0;
    out[432] = 0.0;
    out[433] = 0.0;
    out[434] = 0.0;
    out[435] = 0.0;
    out[436] = 0.0;
    out[437] = 0.0;
    out[438] = 0.0;
    out[439] = 0.0;
    out[440] = 0.0;
    out[441] = 0.0;
    out[442] = 0.0;
    out[443] = 0.0;
    out[444] = 0.0;
    out[445] = 0.0;
    out[446] = 0.0;
    out[447] = 0.0;
    out[448] = 0.0;
    out[449] = 0.0;
    out[450] = 0.0;
    out[451] = 0.0;
    out[452] = 0.0;
    out[453] = 0.0;
    out[454] = 0.0;
    out[455] = 0.0;
    out[456] = 0.0;
    out[457] = 0.0;
    out[458] = 0.0;
    out[459] = 0.0;
    out[460] = 0.0;
    out[461] = 0.0;
    out[462] = 0.0;
    out[463] = 0.0;
    out[464] = 0.0;
    out[465] = 0.0;
    out[466] = 0.0;
    out[467] = 0.0;
    out[468] = 0.0;
    out[469] = 0.0;
    out[470] = 0.0;
    out[471] = 0.0;
    out[472] = 0.0;
    out[473] = 0.0;
    out[474] = 0.0;
    out[475] = 0.0;
    out[476] = 0.0;
    out[477] = 0.0;
    out[478] = 0.0;
    out[479] = 0.0;
    out[480] = 0.0;
    out[481] = 0.0;
    out[482] = 0.0;
    out[483] = 0.0;
    out[484] = 0.0;
    out[485] = 0.0;
    out[486] = 0.0;
    out[487] = 0.0;
    out[488] = 0.0;
    out[489] = t97;
    out[490] = t102;
    out[491] = t104;
    out[492] = 0.0;
    out[493] = 0.0;
    out[494] = 0.0;
    out[495] = 0.0;
    out[496] = 0.0;
    out[497] = 0.0;
    out[498] = t102;
    out[499] = t100;
    out[500] = t88;
    out[501] = 0.0;
    out[502] = 0.0;
    out[503] = 0.0;
    out[504] = 0.0;
    out[505] = 0.0;
    out[506] = 0.0;
    out[507] = t104;
    out[508] = t88;
    out[509] = t99;
    out[510] = 0.0;
    out[511] = 0.0;
    out[512] = 0.0;
    out[513] = 0.0;
    out[514] = 0.0;
    out[515] = 0.0;
    out[516] = 0.0;
    out[517] = 0.0;
    out[518] = 0.0;
    out[519] = 0.0;
    out[520] = 0.0;
    out[521] = 0.0;
    out[522] = 0.0;
    out[523] = 0.0;
    out[524] = 0.0;
    out[525] = 0.0;
    out[526] = 0.0;
    out[527] = 0.0;
    out[528] = 0.0;
    out[529] = 0.0;
    out[530] = 0.0;
    out[531] = 0.0;
    out[532] = 0.0;
    out[533] = 0.0;
    out[534] = 0.0;
    out[535] = 0.0;
    out[536] = 0.0;
    out[537] = 0.0;
    out[538] = 0.0;
    out[539] = 0.0;
    out[540] = 0.0;
    out[541] = 0.0;
    out[542] = 0.0;
    out[543] = -t97;
    out[544] = t108;
    out[545] = t110;
    out[546] = 0.0;
    out[547] = 0.0;
    out[548] = 0.0;
    out[549] = 0.0;
    out[550] = 0.0;
    out[551] = 0.0;
    out[552] = t108;
    out[553] = t106;
    out[554] = t87;
    out[555] = 0.0;
    out[556] = 0.0;
    out[557] = 0.0;
    out[558] = 0.0;
    out[559] = 0.0;
    out[560] = 0.0;
    out[561] = t110;
    out[562] = t87;
    out[563] = t105;
    out[564] = 0.0;
    out[565] = 0.0;
    out[566] = 0.0;
    out[567] = 0.0;
    out[568] = 0.0;
    out[569] = 0.0;
    out[570] = t102;
    out[571] = t100;
    out[572] = t88;
    out[573] = 0.0;
    out[574] = 0.0;
    out[575] = 0.0;
    out[576] = 0.0;
    out[577] = 0.0;
    out[578] = 0.0;
    out[579] = t100;
    out[580] = t96;
    out[581] = t103;
    out[582] = 0.0;
    out[583] = 0.0;
    out[584] = 0.0;
    out[585] = 0.0;
    out[586] = 0.0;
    out[587] = 0.0;
    out[588] = t88;
    out[589] = t103;
    out[590] = t101;
    out[591] = 0.0;
    out[592] = 0.0;
    out[593] = 0.0;
    out[594] = 0.0;
    out[595] = 0.0;
    out[596] = 0.0;
    out[597] = 0.0;
    out[598] = 0.0;
    out[599] = 0.0;
    out[600] = 0.0;
    out[601] = 0.0;
    out[602] = 0.0;
    out[603] = 0.0;
    out[604] = 0.0;
    out[605] = 0.0;
    out[606] = 0.0;
    out[607] = 0.0;
    out[608] = 0.0;
    out[609] = 0.0;
    out[610] = 0.0;
    out[611] = 0.0;
    out[612] = 0.0;
    out[613] = 0.0;
    out[614] = 0.0;
    out[615] = 0.0;
    out[616] = 0.0;
    out[617] = 0.0;
    out[618] = 0.0;
    out[619] = 0.0;
    out[620] = 0.0;
    out[621] = 0.0;
    out[622] = 0.0;
    out[623] = 0.0;
    out[624] = t108;
    out[625] = t106;
    out[626] = t87;
    out[627] = 0.0;
    out[628] = 0.0;
    out[629] = 0.0;
    out[630] = 0.0;
    out[631] = 0.0;
    out[632] = 0.0;
    out[633] = t106;
    out[634] = -t96;
    out[635] = t109;
    out[636] = 0.0;
    out[637] = 0.0;
    out[638] = 0.0;
    out[639] = 0.0;
    out[640] = 0.0;
    out[641] = 0.0;
    out[642] = t87;
    out[643] = t109;
    out[644] = t107;
    out[645] = 0.0;
    out[646] = 0.0;
    out[647] = 0.0;
    out[648] = 0.0;
    out[649] = 0.0;
    out[650] = 0.0;
    out[651] = t104;
    out[652] = t88;
    out[653] = t99;
    out[654] = 0.0;
    out[655] = 0.0;
    out[656] = 0.0;
    out[657] = 0.0;
    out[658] = 0.0;
    out[659] = 0.0;
    out[660] = t88;
    out[661] = t103;
    out[662] = t101;
    out[663] = 0.0;
    out[664] = 0.0;
    out[665] = 0.0;
    out[666] = 0.0;
    out[667] = 0.0;
    out[668] = 0.0;
    out[669] = t99;
    out[670] = t101;
    out[671] = t95;
    out[672] = 0.0;
    out[673] = 0.0;
    out[674] = 0.0;
    out[675] = 0.0;
    out[676] = 0.0;
    out[677] = 0.0;
    out[678] = 0.0;
    out[679] = 0.0;
    out[680] = 0.0;
    out[681] = 0.0;
    out[682] = 0.0;
    out[683] = 0.0;
    out[684] = 0.0;
    out[685] = 0.0;
    out[686] = 0.0;
    out[687] = 0.0;
    out[688] = 0.0;
    out[689] = 0.0;
    out[690] = 0.0;
    out[691] = 0.0;
    out[692] = 0.0;
    out[693] = 0.0;
    out[694] = 0.0;
    out[695] = 0.0;
    out[696] = 0.0;
    out[697] = 0.0;
    out[698] = 0.0;
    out[699] = 0.0;
    out[700] = 0.0;
    out[701] = 0.0;
    out[702] = 0.0;
    out[703] = 0.0;
    out[704] = 0.0;
    out[705] = t110;
    out[706] = t87;
    out[707] = t105;
    out[708] = 0.0;
    out[709] = 0.0;
    out[710] = 0.0;
    out[711] = 0.0;
    out[712] = 0.0;
    out[713] = 0.0;
    out[714] = t87;
    out[715] = t109;
    out[716] = t107;
    out[717] = 0.0;
    out[718] = 0.0;
    out[719] = 0.0;
    out[720] = 0.0;
    out[721] = 0.0;
    out[722] = 0.0;
    out[723] = t105;
    out[724] = t107;
    out[725] = -t95;
    out[726] = 0.0;
    out[727] = 0.0;
    out[728] = 0.0;
}


/** Jacobian d(df/dx)/du in vector form (column-wise **/
void Chain_2_ProblemDescription::dfdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
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
    out[43] = 0.0;
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
    out[64] = 0.0;
    out[65] = 0.0;
    out[66] = 0.0;
    out[67] = 0.0;
    out[68] = 0.0;
    out[69] = 0.0;
    out[70] = 0.0;
    out[71] = 0.0;
    out[72] = 0.0;
    out[73] = 0.0;
    out[74] = 0.0;
    out[75] = 0.0;
    out[76] = 0.0;
    out[77] = 0.0;
    out[78] = 0.0;
    out[79] = 0.0;
    out[80] = 0.0;
    out[81] = 0.0;
    out[82] = 0.0;
    out[83] = 0.0;
    out[84] = 0.0;
    out[85] = 0.0;
    out[86] = 0.0;
    out[87] = 0.0;
    out[88] = 0.0;
    out[89] = 0.0;
    out[90] = 0.0;
    out[91] = 0.0;
    out[92] = 0.0;
    out[93] = 0.0;
    out[94] = 0.0;
    out[95] = 0.0;
    out[96] = 0.0;
    out[97] = 0.0;
    out[98] = 0.0;
    out[99] = 0.0;
    out[100] = 0.0;
    out[101] = 0.0;
    out[102] = 0.0;
    out[103] = 0.0;
    out[104] = 0.0;
    out[105] = 0.0;
    out[106] = 0.0;
    out[107] = 0.0;
    out[108] = 0.0;
    out[109] = 0.0;
    out[110] = 0.0;
    out[111] = 0.0;
    out[112] = 0.0;
    out[113] = 0.0;
    out[114] = 0.0;
    out[115] = 0.0;
    out[116] = 0.0;
    out[117] = 0.0;
    out[118] = 0.0;
    out[119] = 0.0;
    out[120] = 0.0;
    out[121] = 0.0;
    out[122] = 0.0;
    out[123] = 0.0;
    out[124] = 0.0;
    out[125] = 0.0;
    out[126] = 0.0;
    out[127] = 0.0;
    out[128] = 0.0;
    out[129] = 0.0;
    out[130] = 0.0;
    out[131] = 0.0;
    out[132] = 0.0;
    out[133] = 0.0;
    out[134] = 0.0;
    out[135] = 0.0;
    out[136] = 0.0;
    out[137] = 0.0;
    out[138] = 0.0;
    out[139] = 0.0;
    out[140] = 0.0;
    out[141] = 0.0;
    out[142] = 0.0;
    out[143] = 0.0;
    out[144] = 0.0;
    out[145] = 0.0;
    out[146] = 0.0;
    out[147] = 0.0;
    out[148] = 0.0;
    out[149] = 0.0;
    out[150] = 0.0;
    out[151] = 0.0;
    out[152] = 0.0;
    out[153] = 0.0;
    out[154] = 0.0;
    out[155] = 0.0;
    out[156] = 0.0;
    out[157] = 0.0;
    out[158] = 0.0;
    out[159] = 0.0;
    out[160] = 0.0;
    out[161] = 0.0;
    out[162] = 0.0;
    out[163] = 0.0;
    out[164] = 0.0;
    out[165] = 0.0;
    out[166] = 0.0;
    out[167] = 0.0;
    out[168] = 0.0;
    out[169] = 0.0;
    out[170] = 0.0;
    out[171] = 0.0;
    out[172] = 0.0;
    out[173] = 0.0;
    out[174] = 0.0;
    out[175] = 0.0;
    out[176] = 0.0;
    out[177] = 0.0;
    out[178] = 0.0;
    out[179] = 0.0;
    out[180] = 0.0;
    out[181] = 0.0;
    out[182] = 0.0;
    out[183] = 0.0;
    out[184] = 0.0;
    out[185] = 0.0;
    out[186] = 0.0;
    out[187] = 0.0;
    out[188] = 0.0;
    out[189] = 0.0;
    out[190] = 0.0;
    out[191] = 0.0;
    out[192] = 0.0;
    out[193] = 0.0;
    out[194] = 0.0;
    out[195] = 0.0;
    out[196] = 0.0;
    out[197] = 0.0;
    out[198] = 0.0;
    out[199] = 0.0;
    out[200] = 0.0;
    out[201] = 0.0;
    out[202] = 0.0;
    out[203] = 0.0;
    out[204] = 0.0;
    out[205] = 0.0;
    out[206] = 0.0;
    out[207] = 0.0;
    out[208] = 0.0;
    out[209] = 0.0;
    out[210] = 0.0;
    out[211] = 0.0;
    out[212] = 0.0;
    out[213] = 0.0;
    out[214] = 0.0;
    out[215] = 0.0;
    out[216] = 0.0;
    out[217] = 0.0;
    out[218] = 0.0;
    out[219] = 0.0;
    out[220] = 0.0;
    out[221] = 0.0;
    out[222] = 0.0;
    out[223] = 0.0;
    out[224] = 0.0;
    out[225] = 0.0;
    out[226] = 0.0;
    out[227] = 0.0;
    out[228] = 0.0;
    out[229] = 0.0;
    out[230] = 0.0;
    out[231] = 0.0;
    out[232] = 0.0;
    out[233] = 0.0;
    out[234] = 0.0;
    out[235] = 0.0;
    out[236] = 0.0;
    out[237] = 0.0;
    out[238] = 0.0;
    out[239] = 0.0;
    out[240] = 0.0;
    out[241] = 0.0;
    out[242] = 0.0;
}


/** Jacobian d(df/dx)/dp in vector form (column-wise **/
void Chain_2_ProblemDescription::dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian df/dp in vector form (column-wise **/
void Chain_2_ProblemDescription::dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian d(df/dp)/du in vector form (column-wise **/
void Chain_2_ProblemDescription::dfdpdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}

/** Jacobian dh/dx in vector form (column-wise) **/
void Chain_2_ProblemDescription::dhdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian dh/du in vector form (column-wise) **/
void Chain_2_ProblemDescription::dhdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Hessian d(dh/dx)/dx in vector form (column-wise) **/
void Chain_2_ProblemDescription::dhdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian d(dh/dx)/du in vector form (column-wise) **/
void Chain_2_ProblemDescription::dhdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}