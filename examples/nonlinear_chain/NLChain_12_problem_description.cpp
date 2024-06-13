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
* This probfct-file describes the nonlinear chain problem with 12 chain elements from
* Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
* Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
* Process Control 22(3), 540-550 (2012)
*
*/

#include "NLChain_12_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

using namespace grampc;


Chain_12_ProblemDescription::Chain_12_ProblemDescription(const std::vector<typeRNum>& pCost)
: pCost_(pCost)
{
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void Chain_12_ProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 69;
    *Nu  = 3;
    *Np  = 0;
    *Nh  = 0;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void Chain_12_ProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
    ctypeRNum  t2 = x[0]*x[0];
    ctypeRNum  t3 = x[1]*x[1];
    ctypeRNum  t4 = x[2]*x[2];
    ctypeRNum  t5 = -x[6];
    ctypeRNum  t6 = -x[7];
    ctypeRNum  t7 = -x[8];
    ctypeRNum  t8 = -x[12];
    ctypeRNum  t9 = -x[13];
    ctypeRNum  t10 = -x[14];
    ctypeRNum  t11 = -x[18];
    ctypeRNum  t12 = -x[19];
    ctypeRNum  t13 = -x[20];
    ctypeRNum  t14 = -x[24];
    ctypeRNum  t15 = -x[25];
    ctypeRNum  t16 = -x[26];
    ctypeRNum  t17 = -x[30];
    ctypeRNum  t18 = -x[31];
    ctypeRNum  t19 = -x[32];
    ctypeRNum  t20 = -x[36];
    ctypeRNum  t21 = -x[37];
    ctypeRNum  t22 = -x[38];
    ctypeRNum  t23 = -x[42];
    ctypeRNum  t24 = -x[43];
    ctypeRNum  t25 = -x[44];
    ctypeRNum  t26 = -x[48];
    ctypeRNum  t27 = -x[49];
    ctypeRNum  t28 = -x[50];
    ctypeRNum  t29 = -x[54];
    ctypeRNum  t30 = -x[55];
    ctypeRNum  t31 = -x[56];
    ctypeRNum  t32 = -x[60];
    ctypeRNum  t33 = -x[61];
    ctypeRNum  t34 = -x[62];
    ctypeRNum  t35 = -x[66];
    ctypeRNum  t36 = -x[67];
    ctypeRNum  t37 = -x[68];
    ctypeRNum  t38 = t5+x[0];
    ctypeRNum  t39 = t6+x[1];
    ctypeRNum  t40 = t7+x[2];
    ctypeRNum  t41 = t8+x[6];
    ctypeRNum  t42 = t9+x[7];
    ctypeRNum  t43 = t10+x[8];
    ctypeRNum  t44 = t11+x[12];
    ctypeRNum  t45 = t12+x[13];
    ctypeRNum  t46 = t13+x[14];
    ctypeRNum  t47 = t14+x[18];
    ctypeRNum  t48 = t15+x[19];
    ctypeRNum  t49 = t16+x[20];
    ctypeRNum  t50 = t17+x[24];
    ctypeRNum  t51 = t18+x[25];
    ctypeRNum  t52 = t19+x[26];
    ctypeRNum  t53 = t20+x[30];
    ctypeRNum  t54 = t21+x[31];
    ctypeRNum  t55 = t22+x[32];
    ctypeRNum  t56 = t23+x[36];
    ctypeRNum  t57 = t24+x[37];
    ctypeRNum  t58 = t25+x[38];
    ctypeRNum  t59 = t26+x[42];
    ctypeRNum  t60 = t27+x[43];
    ctypeRNum  t61 = t28+x[44];
    ctypeRNum  t62 = t29+x[48];
    ctypeRNum  t63 = t30+x[49];
    ctypeRNum  t64 = t31+x[50];
    ctypeRNum  t65 = t32+x[54];
    ctypeRNum  t66 = t33+x[55];
    ctypeRNum  t67 = t34+x[56];
    ctypeRNum  t68 = t35+x[60];
    ctypeRNum  t69 = t36+x[61];
    ctypeRNum  t70 = t37+x[62];
    ctypeRNum  t104 = t2+t3+t4;
    ctypeRNum  t71 = t38*t38;
    ctypeRNum  t72 = t39*t39;
    ctypeRNum  t73 = t40*t40;
    ctypeRNum  t74 = t41*t41;
    ctypeRNum  t75 = t42*t42;
    ctypeRNum  t76 = t43*t43;
    ctypeRNum  t77 = t44*t44;
    ctypeRNum  t78 = t45*t45;
    ctypeRNum  t79 = t46*t46;
    ctypeRNum  t80 = t47*t47;
    ctypeRNum  t81 = t48*t48;
    ctypeRNum  t82 = t49*t49;
    ctypeRNum  t83 = t50*t50;
    ctypeRNum  t84 = t51*t51;
    ctypeRNum  t85 = t52*t52;
    ctypeRNum  t86 = t53*t53;
    ctypeRNum  t87 = t54*t54;
    ctypeRNum  t88 = t55*t55;
    ctypeRNum  t89 = t56*t56;
    ctypeRNum  t90 = t57*t57;
    ctypeRNum  t91 = t58*t58;
    ctypeRNum  t92 = t59*t59;
    ctypeRNum  t93 = t60*t60;
    ctypeRNum  t94 = t61*t61;
    ctypeRNum  t95 = t62*t62;
    ctypeRNum  t96 = t63*t63;
    ctypeRNum  t97 = t64*t64;
    ctypeRNum  t98 = t65*t65;
    ctypeRNum  t99 = t66*t66;
    ctypeRNum  t100 = t67*t67;
    ctypeRNum  t101 = t68*t68;
    ctypeRNum  t102 = t69*t69;
    ctypeRNum  t103 = t70*t70;
    ctypeRNum  t105 = 1.0/sqrt(t104);
    ctypeRNum  t106 = t105*(1.1E+1/4.0E+2);
    ctypeRNum  t107 = t71+t72+t73;
    ctypeRNum  t108 = t74+t75+t76;
    ctypeRNum  t109 = t77+t78+t79;
    ctypeRNum  t110 = t80+t81+t82;
    ctypeRNum  t111 = t83+t84+t85;
    ctypeRNum  t112 = t86+t87+t88;
    ctypeRNum  t113 = t89+t90+t91;
    ctypeRNum  t114 = t92+t93+t94;
    ctypeRNum  t115 = t95+t96+t97;
    ctypeRNum  t116 = t98+t99+t100;
    ctypeRNum  t117 = t101+t102+t103;
    ctypeRNum  t118 = t106-3.0/5.0;
    ctypeRNum  t119 = 1.0/sqrt(t107);
    ctypeRNum  t120 = 1.0/sqrt(t108);
    ctypeRNum  t121 = 1.0/sqrt(t109);
    ctypeRNum  t122 = 1.0/sqrt(t110);
    ctypeRNum  t123 = 1.0/sqrt(t111);
    ctypeRNum  t124 = 1.0/sqrt(t112);
    ctypeRNum  t125 = 1.0/sqrt(t113);
    ctypeRNum  t126 = 1.0/sqrt(t114);
    ctypeRNum  t127 = 1.0/sqrt(t115);
    ctypeRNum  t128 = 1.0/sqrt(t116);
    ctypeRNum  t129 = 1.0/sqrt(t117);
    ctypeRNum  t130 = t119*(1.1E+1/4.0E+2);
    ctypeRNum  t131 = t120*(1.1E+1/4.0E+2);
    ctypeRNum  t132 = t121*(1.1E+1/4.0E+2);
    ctypeRNum  t133 = t122*(1.1E+1/4.0E+2);
    ctypeRNum  t134 = t123*(1.1E+1/4.0E+2);
    ctypeRNum  t135 = t124*(1.1E+1/4.0E+2);
    ctypeRNum  t136 = t125*(1.1E+1/4.0E+2);
    ctypeRNum  t137 = t126*(1.1E+1/4.0E+2);
    ctypeRNum  t138 = t127*(1.1E+1/4.0E+2);
    ctypeRNum  t139 = t128*(1.1E+1/4.0E+2);
    ctypeRNum  t140 = t129*(1.1E+1/4.0E+2);
    ctypeRNum  t141 = t130-3.0/5.0;
    ctypeRNum  t142 = t131-3.0/5.0;
    ctypeRNum  t143 = t132-3.0/5.0;
    ctypeRNum  t144 = t133-3.0/5.0;
    ctypeRNum  t145 = t134-3.0/5.0;
    ctypeRNum  t146 = t135-3.0/5.0;
    ctypeRNum  t147 = t136-3.0/5.0;
    ctypeRNum  t148 = t137-3.0/5.0;
    ctypeRNum  t149 = t138-3.0/5.0;
    ctypeRNum  t150 = t139-3.0/5.0;
    ctypeRNum  t151 = t140-3.0/5.0;
    ctypeRNum  t152 = t38*t141*(1.0E+2/3.0);
    ctypeRNum  t153 = t39*t141*(1.0E+2/3.0);
    ctypeRNum  t154 = t40*t141*(1.0E+2/3.0);
    ctypeRNum  t155 = t41*t142*(1.0E+2/3.0);
    ctypeRNum  t156 = t42*t142*(1.0E+2/3.0);
    ctypeRNum  t157 = t43*t142*(1.0E+2/3.0);
    ctypeRNum  t158 = t44*t143*(1.0E+2/3.0);
    ctypeRNum  t159 = t45*t143*(1.0E+2/3.0);
    ctypeRNum  t160 = t46*t143*(1.0E+2/3.0);
    ctypeRNum  t161 = t47*t144*(1.0E+2/3.0);
    ctypeRNum  t162 = t48*t144*(1.0E+2/3.0);
    ctypeRNum  t163 = t49*t144*(1.0E+2/3.0);
    ctypeRNum  t164 = t50*t145*(1.0E+2/3.0);
    ctypeRNum  t165 = t51*t145*(1.0E+2/3.0);
    ctypeRNum  t166 = t52*t145*(1.0E+2/3.0);
    ctypeRNum  t167 = t53*t146*(1.0E+2/3.0);
    ctypeRNum  t168 = t54*t146*(1.0E+2/3.0);
    ctypeRNum  t169 = t55*t146*(1.0E+2/3.0);
    ctypeRNum  t170 = t56*t147*(1.0E+2/3.0);
    ctypeRNum  t171 = t57*t147*(1.0E+2/3.0);
    ctypeRNum  t172 = t58*t147*(1.0E+2/3.0);
    ctypeRNum  t173 = t59*t148*(1.0E+2/3.0);
    ctypeRNum  t174 = t60*t148*(1.0E+2/3.0);
    ctypeRNum  t175 = t61*t148*(1.0E+2/3.0);
    ctypeRNum  t176 = t62*t149*(1.0E+2/3.0);
    ctypeRNum  t177 = t63*t149*(1.0E+2/3.0);
    ctypeRNum  t178 = t64*t149*(1.0E+2/3.0);
    ctypeRNum  t179 = t65*t150*(1.0E+2/3.0);
    ctypeRNum  t180 = t66*t150*(1.0E+2/3.0);
    ctypeRNum  t181 = t67*t150*(1.0E+2/3.0);
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = t152+t118*x[0]*(1.0E+2/3.0);
    out[4] = t153+t118*x[1]*(1.0E+2/3.0);
    out[5] = t154+t118*x[2]*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[6] = x[9];
    out[7] = x[10];
    out[8] = x[11];
    out[9] = -t152+t155;
    out[10] = -t153+t156;
    out[11] = -t154+t157-9.81E+2/1.0E+2;
    out[12] = x[15];
    out[13] = x[16];
    out[14] = x[17];
    out[15] = -t155+t158;
    out[16] = -t156+t159;
    out[17] = -t157+t160-9.81E+2/1.0E+2;
    out[18] = x[21];
    out[19] = x[22];
    out[20] = x[23];
    out[21] = -t158+t161;
    out[22] = -t159+t162;
    out[23] = -t160+t163-9.81E+2/1.0E+2;
    out[24] = x[27];
    out[25] = x[28];
    out[26] = x[29];
    out[27] = -t161+t164;
    out[28] = -t162+t165;
    out[29] = -t163+t166-9.81E+2/1.0E+2;
    out[30] = x[33];
    out[31] = x[34];
    out[32] = x[35];
    out[33] = -t164+t167;
    out[34] = -t165+t168;
    out[35] = -t166+t169-9.81E+2/1.0E+2;
    out[36] = x[39];
    out[37] = x[40];
    out[38] = x[41];
    out[39] = -t167+t170;
    out[40] = -t168+t171;
    out[41] = -t169+t172-9.81E+2/1.0E+2;
    out[42] = x[45];
    out[43] = x[46];
    out[44] = x[47];
    out[45] = -t170+t173;
    out[46] = -t171+t174;
    out[47] = -t172+t175-9.81E+2/1.0E+2;
    out[48] = x[51];
    out[49] = x[52];
    out[50] = x[53];
    out[51] = -t173+t176;
    out[52] = -t174+t177;
    out[53] = -t175+t178-9.81E+2/1.0E+2;
    out[54] = x[57];
    out[55] = x[58];
    out[56] = x[59];
    out[57] = -t176+t179;
    out[58] = -t177+t180;
    out[59] = -t178+t181-9.81E+2/1.0E+2;
    out[60] = x[63];
    out[61] = x[64];
    out[62] = x[65];
    out[63] = -t179+t68*t151*(1.0E+2/3.0);
    out[64] = -t180+t69*t151*(1.0E+2/3.0);
    out[65] = -t181+t70*t151*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[66] = u[0];
    out[67] = u[1];
    out[68] = u[2];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Chain_12_ProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
    ctypeRNum  t2 = x[0]*2.0;
    ctypeRNum  t3 = x[1]*2.0;
    ctypeRNum  t4 = x[2]*2.0;
    ctypeRNum  t5 = x[6]*2.0;
    ctypeRNum  t6 = x[7]*2.0;
    ctypeRNum  t7 = x[8]*2.0;
    ctypeRNum  t8 = x[12]*2.0;
    ctypeRNum  t9 = x[13]*2.0;
    ctypeRNum  t10 = x[14]*2.0;
    ctypeRNum  t11 = x[18]*2.0;
    ctypeRNum  t12 = x[19]*2.0;
    ctypeRNum  t13 = x[20]*2.0;
    ctypeRNum  t14 = x[24]*2.0;
    ctypeRNum  t15 = x[25]*2.0;
    ctypeRNum  t16 = x[26]*2.0;
    ctypeRNum  t17 = x[30]*2.0;
    ctypeRNum  t18 = x[31]*2.0;
    ctypeRNum  t19 = x[32]*2.0;
    ctypeRNum  t20 = x[36]*2.0;
    ctypeRNum  t21 = x[37]*2.0;
    ctypeRNum  t22 = x[38]*2.0;
    ctypeRNum  t23 = x[42]*2.0;
    ctypeRNum  t24 = x[43]*2.0;
    ctypeRNum  t25 = x[44]*2.0;
    ctypeRNum  t26 = x[48]*2.0;
    ctypeRNum  t27 = x[49]*2.0;
    ctypeRNum  t28 = x[50]*2.0;
    ctypeRNum  t29 = x[54]*2.0;
    ctypeRNum  t30 = x[55]*2.0;
    ctypeRNum  t31 = x[56]*2.0;
    ctypeRNum  t32 = x[60]*2.0;
    ctypeRNum  t33 = x[61]*2.0;
    ctypeRNum  t34 = x[62]*2.0;
    ctypeRNum  t35 = x[66]*2.0;
    ctypeRNum  t36 = x[67]*2.0;
    ctypeRNum  t37 = x[68]*2.0;
    ctypeRNum  t38 = x[0]*x[0];
    ctypeRNum  t39 = x[1]*x[1];
    ctypeRNum  t40 = x[2]*x[2];
    ctypeRNum  t41 = x[60]*x[60];
    ctypeRNum  t42 = x[61]*x[61];
    ctypeRNum  t43 = x[62]*x[62];
    ctypeRNum  t44 = -x[6];
    ctypeRNum  t46 = -x[7];
    ctypeRNum  t48 = -x[8];
    ctypeRNum  t50 = -x[12];
    ctypeRNum  t52 = -x[13];
    ctypeRNum  t54 = -x[14];
    ctypeRNum  t56 = -x[18];
    ctypeRNum  t58 = -x[19];
    ctypeRNum  t60 = -x[20];
    ctypeRNum  t62 = -x[24];
    ctypeRNum  t64 = -x[25];
    ctypeRNum  t66 = -x[26];
    ctypeRNum  t68 = -x[30];
    ctypeRNum  t70 = -x[31];
    ctypeRNum  t72 = -x[32];
    ctypeRNum  t74 = -x[36];
    ctypeRNum  t76 = -x[37];
    ctypeRNum  t78 = -x[38];
    ctypeRNum  t80 = -x[42];
    ctypeRNum  t82 = -x[43];
    ctypeRNum  t84 = -x[44];
    ctypeRNum  t86 = -x[48];
    ctypeRNum  t88 = -x[49];
    ctypeRNum  t90 = -x[50];
    ctypeRNum  t92 = -x[54];
    ctypeRNum  t94 = -x[55];
    ctypeRNum  t96 = -x[56];
    ctypeRNum  t98 = -x[60];
    ctypeRNum  t100 = -x[61];
    ctypeRNum  t102 = -x[62];
    ctypeRNum  t104 = -x[66];
    ctypeRNum  t106 = -x[67];
    ctypeRNum  t108 = -x[68];
    ctypeRNum  t110 = vec[63]*x[60]*1.1E+1;
    ctypeRNum  t111 = vec[64]*x[61]*1.1E+1;
    ctypeRNum  t112 = vec[65]*x[62]*1.1E+1;
    ctypeRNum  t113 = vec[63]*x[66]*1.1E+1;
    ctypeRNum  t114 = vec[64]*x[67]*1.1E+1;
    ctypeRNum  t115 = vec[65]*x[68]*1.1E+1;
    ctypeRNum  t45 = -t5;
    ctypeRNum  t47 = -t6;
    ctypeRNum  t49 = -t7;
    ctypeRNum  t51 = -t8;
    ctypeRNum  t53 = -t9;
    ctypeRNum  t55 = -t10;
    ctypeRNum  t57 = -t11;
    ctypeRNum  t59 = -t12;
    ctypeRNum  t61 = -t13;
    ctypeRNum  t63 = -t14;
    ctypeRNum  t65 = -t15;
    ctypeRNum  t67 = -t16;
    ctypeRNum  t69 = -t17;
    ctypeRNum  t71 = -t18;
    ctypeRNum  t73 = -t19;
    ctypeRNum  t75 = -t20;
    ctypeRNum  t77 = -t21;
    ctypeRNum  t79 = -t22;
    ctypeRNum  t81 = -t23;
    ctypeRNum  t83 = -t24;
    ctypeRNum  t85 = -t25;
    ctypeRNum  t87 = -t26;
    ctypeRNum  t89 = -t27;
    ctypeRNum  t91 = -t28;
    ctypeRNum  t93 = -t29;
    ctypeRNum  t95 = -t30;
    ctypeRNum  t97 = -t31;
    ctypeRNum  t99 = -t32;
    ctypeRNum  t101 = -t33;
    ctypeRNum  t103 = -t34;
    ctypeRNum  t105 = -t35;
    ctypeRNum  t107 = -t36;
    ctypeRNum  t109 = -t37;
    ctypeRNum  t116 = -t113;
    ctypeRNum  t117 = -t114;
    ctypeRNum  t118 = -t115;
    ctypeRNum  t119 = t44+x[0];
    ctypeRNum  t120 = t46+x[1];
    ctypeRNum  t121 = t48+x[2];
    ctypeRNum  t122 = t50+x[6];
    ctypeRNum  t123 = t52+x[7];
    ctypeRNum  t124 = t54+x[8];
    ctypeRNum  t125 = t56+x[12];
    ctypeRNum  t126 = t58+x[13];
    ctypeRNum  t127 = t60+x[14];
    ctypeRNum  t128 = t62+x[18];
    ctypeRNum  t129 = t64+x[19];
    ctypeRNum  t130 = t66+x[20];
    ctypeRNum  t131 = t68+x[24];
    ctypeRNum  t132 = t70+x[25];
    ctypeRNum  t133 = t72+x[26];
    ctypeRNum  t134 = t74+x[30];
    ctypeRNum  t135 = t76+x[31];
    ctypeRNum  t136 = t78+x[32];
    ctypeRNum  t137 = t80+x[36];
    ctypeRNum  t138 = t82+x[37];
    ctypeRNum  t139 = t84+x[38];
    ctypeRNum  t140 = t86+x[42];
    ctypeRNum  t141 = t88+x[43];
    ctypeRNum  t142 = t90+x[44];
    ctypeRNum  t143 = t92+x[48];
    ctypeRNum  t144 = t94+x[49];
    ctypeRNum  t145 = t96+x[50];
    ctypeRNum  t146 = t98+x[54];
    ctypeRNum  t147 = t100+x[55];
    ctypeRNum  t148 = t102+x[56];
    ctypeRNum  t149 = t104+x[60];
    ctypeRNum  t150 = t106+x[61];
    ctypeRNum  t151 = t108+x[62];
    ctypeRNum  t218 = t38+t39+t40;
    ctypeRNum  t152 = t2+t45;
    ctypeRNum  t153 = t3+t47;
    ctypeRNum  t154 = t4+t49;
    ctypeRNum  t155 = t5+t51;
    ctypeRNum  t156 = t6+t53;
    ctypeRNum  t157 = t7+t55;
    ctypeRNum  t158 = t8+t57;
    ctypeRNum  t159 = t9+t59;
    ctypeRNum  t160 = t10+t61;
    ctypeRNum  t161 = t11+t63;
    ctypeRNum  t162 = t12+t65;
    ctypeRNum  t163 = t13+t67;
    ctypeRNum  t164 = t14+t69;
    ctypeRNum  t165 = t15+t71;
    ctypeRNum  t166 = t16+t73;
    ctypeRNum  t167 = t17+t75;
    ctypeRNum  t168 = t18+t77;
    ctypeRNum  t169 = t19+t79;
    ctypeRNum  t170 = t20+t81;
    ctypeRNum  t171 = t21+t83;
    ctypeRNum  t172 = t22+t85;
    ctypeRNum  t173 = t23+t87;
    ctypeRNum  t174 = t24+t89;
    ctypeRNum  t175 = t25+t91;
    ctypeRNum  t176 = t26+t93;
    ctypeRNum  t177 = t27+t95;
    ctypeRNum  t178 = t28+t97;
    ctypeRNum  t179 = t29+t99;
    ctypeRNum  t180 = t30+t101;
    ctypeRNum  t181 = t31+t103;
    ctypeRNum  t182 = t32+t105;
    ctypeRNum  t183 = t33+t107;
    ctypeRNum  t184 = t34+t109;
    ctypeRNum  t185 = t119*t119;
    ctypeRNum  t186 = t120*t120;
    ctypeRNum  t187 = t121*t121;
    ctypeRNum  t188 = t122*t122;
    ctypeRNum  t189 = t123*t123;
    ctypeRNum  t190 = t124*t124;
    ctypeRNum  t191 = t125*t125;
    ctypeRNum  t192 = t126*t126;
    ctypeRNum  t193 = t127*t127;
    ctypeRNum  t194 = t128*t128;
    ctypeRNum  t195 = t129*t129;
    ctypeRNum  t196 = t130*t130;
    ctypeRNum  t197 = t131*t131;
    ctypeRNum  t198 = t132*t132;
    ctypeRNum  t199 = t133*t133;
    ctypeRNum  t200 = t134*t134;
    ctypeRNum  t201 = t135*t135;
    ctypeRNum  t202 = t136*t136;
    ctypeRNum  t203 = t137*t137;
    ctypeRNum  t204 = t138*t138;
    ctypeRNum  t205 = t139*t139;
    ctypeRNum  t206 = t140*t140;
    ctypeRNum  t207 = t141*t141;
    ctypeRNum  t208 = t142*t142;
    ctypeRNum  t209 = t143*t143;
    ctypeRNum  t210 = t144*t144;
    ctypeRNum  t211 = t145*t145;
    ctypeRNum  t212 = t146*t146;
    ctypeRNum  t213 = t147*t147;
    ctypeRNum  t214 = t148*t148;
    ctypeRNum  t215 = t149*t149;
    ctypeRNum  t216 = t150*t150;
    ctypeRNum  t217 = t151*t151;
    ctypeRNum  t219 = 1.0/sqrt(t218);
    ctypeRNum  t220 = t219*t219*t219;
    ctypeRNum  t221 = t219*(1.1E+1/1.2E+1);
    ctypeRNum  t226 = t185+t186+t187;
    ctypeRNum  t227 = t188+t189+t190;
    ctypeRNum  t228 = t191+t192+t193;
    ctypeRNum  t229 = t194+t195+t196;
    ctypeRNum  t230 = t197+t198+t199;
    ctypeRNum  t231 = t200+t201+t202;
    ctypeRNum  t232 = t203+t204+t205;
    ctypeRNum  t233 = t206+t207+t208;
    ctypeRNum  t234 = t209+t210+t211;
    ctypeRNum  t235 = t212+t213+t214;
    ctypeRNum  t236 = t215+t216+t217;
    ctypeRNum  t222 = -t221;
    ctypeRNum  t223 = t220*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t224 = t220*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t225 = t220*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t237 = pow(t236,3.0/2.0);
    ctypeRNum  t238 = 1.0/sqrt(t226);
    ctypeRNum  t240 = 1.0/sqrt(t227);
    ctypeRNum  t242 = 1.0/sqrt(t228);
    ctypeRNum  t244 = 1.0/sqrt(t229);
    ctypeRNum  t246 = 1.0/sqrt(t230);
    ctypeRNum  t248 = 1.0/sqrt(t231);
    ctypeRNum  t250 = 1.0/sqrt(t232);
    ctypeRNum  t252 = 1.0/sqrt(t233);
    ctypeRNum  t254 = 1.0/sqrt(t234);
    ctypeRNum  t256 = 1.0/sqrt(t235);
    ctypeRNum  t258 = 1.0/sqrt(t236);
    ctypeRNum  t239 = t238*t238*t238;
    ctypeRNum  t241 = t240*t240*t240;
    ctypeRNum  t243 = t242*t242*t242;
    ctypeRNum  t245 = t244*t244*t244;
    ctypeRNum  t247 = t246*t246*t246;
    ctypeRNum  t249 = t248*t248*t248;
    ctypeRNum  t251 = t250*t250*t250;
    ctypeRNum  t253 = t252*t252*t252;
    ctypeRNum  t255 = t254*t254*t254;
    ctypeRNum  t257 = t256*t256*t256;
    ctypeRNum  t259 = 1.0/t237;
    ctypeRNum  t260 = t238*(1.1E+1/1.2E+1);
    ctypeRNum  t261 = t240*(1.1E+1/1.2E+1);
    ctypeRNum  t262 = t242*(1.1E+1/1.2E+1);
    ctypeRNum  t263 = t244*(1.1E+1/1.2E+1);
    ctypeRNum  t264 = t246*(1.1E+1/1.2E+1);
    ctypeRNum  t265 = t248*(1.1E+1/1.2E+1);
    ctypeRNum  t266 = t250*(1.1E+1/1.2E+1);
    ctypeRNum  t267 = t252*(1.1E+1/1.2E+1);
    ctypeRNum  t268 = t254*(1.1E+1/1.2E+1);
    ctypeRNum  t269 = t256*(1.1E+1/1.2E+1);
    ctypeRNum  t270 = -t260;
    ctypeRNum  t271 = -t261;
    ctypeRNum  t272 = -t262;
    ctypeRNum  t273 = -t263;
    ctypeRNum  t274 = -t264;
    ctypeRNum  t275 = -t265;
    ctypeRNum  t276 = -t266;
    ctypeRNum  t277 = -t267;
    ctypeRNum  t278 = -t268;
    ctypeRNum  t279 = -t269;
    ctypeRNum  t280 = t119*t120*t239*(1.1E+1/1.2E+1);
    ctypeRNum  t281 = t119*t121*t239*(1.1E+1/1.2E+1);
    ctypeRNum  t282 = t120*t121*t239*(1.1E+1/1.2E+1);
    ctypeRNum  t283 = t122*t123*t241*(1.1E+1/1.2E+1);
    ctypeRNum  t284 = t122*t124*t241*(1.1E+1/1.2E+1);
    ctypeRNum  t285 = t123*t124*t241*(1.1E+1/1.2E+1);
    ctypeRNum  t286 = t125*t126*t243*(1.1E+1/1.2E+1);
    ctypeRNum  t287 = t125*t127*t243*(1.1E+1/1.2E+1);
    ctypeRNum  t288 = t126*t127*t243*(1.1E+1/1.2E+1);
    ctypeRNum  t289 = t128*t129*t245*(1.1E+1/1.2E+1);
    ctypeRNum  t290 = t128*t130*t245*(1.1E+1/1.2E+1);
    ctypeRNum  t291 = t129*t130*t245*(1.1E+1/1.2E+1);
    ctypeRNum  t292 = t131*t132*t247*(1.1E+1/1.2E+1);
    ctypeRNum  t293 = t131*t133*t247*(1.1E+1/1.2E+1);
    ctypeRNum  t294 = t132*t133*t247*(1.1E+1/1.2E+1);
    ctypeRNum  t295 = t134*t135*t249*(1.1E+1/1.2E+1);
    ctypeRNum  t296 = t134*t136*t249*(1.1E+1/1.2E+1);
    ctypeRNum  t297 = t135*t136*t249*(1.1E+1/1.2E+1);
    ctypeRNum  t298 = t137*t138*t251*(1.1E+1/1.2E+1);
    ctypeRNum  t299 = t137*t139*t251*(1.1E+1/1.2E+1);
    ctypeRNum  t300 = t138*t139*t251*(1.1E+1/1.2E+1);
    ctypeRNum  t301 = t140*t141*t253*(1.1E+1/1.2E+1);
    ctypeRNum  t302 = t140*t142*t253*(1.1E+1/1.2E+1);
    ctypeRNum  t303 = t141*t142*t253*(1.1E+1/1.2E+1);
    ctypeRNum  t304 = t143*t144*t255*(1.1E+1/1.2E+1);
    ctypeRNum  t305 = t143*t145*t255*(1.1E+1/1.2E+1);
    ctypeRNum  t306 = t144*t145*t255*(1.1E+1/1.2E+1);
    ctypeRNum  t307 = t146*t147*t257*(1.1E+1/1.2E+1);
    ctypeRNum  t308 = t146*t148*t257*(1.1E+1/1.2E+1);
    ctypeRNum  t309 = t147*t148*t257*(1.1E+1/1.2E+1);
    ctypeRNum  t310 = t119*t152*t239*(1.1E+1/2.4E+1);
    ctypeRNum  t311 = t120*t153*t239*(1.1E+1/2.4E+1);
    ctypeRNum  t312 = t121*t154*t239*(1.1E+1/2.4E+1);
    ctypeRNum  t313 = t122*t155*t241*(1.1E+1/2.4E+1);
    ctypeRNum  t314 = t123*t156*t241*(1.1E+1/2.4E+1);
    ctypeRNum  t315 = t124*t157*t241*(1.1E+1/2.4E+1);
    ctypeRNum  t316 = t125*t158*t243*(1.1E+1/2.4E+1);
    ctypeRNum  t317 = t126*t159*t243*(1.1E+1/2.4E+1);
    ctypeRNum  t318 = t127*t160*t243*(1.1E+1/2.4E+1);
    ctypeRNum  t319 = t128*t161*t245*(1.1E+1/2.4E+1);
    ctypeRNum  t320 = t129*t162*t245*(1.1E+1/2.4E+1);
    ctypeRNum  t321 = t130*t163*t245*(1.1E+1/2.4E+1);
    ctypeRNum  t322 = t131*t164*t247*(1.1E+1/2.4E+1);
    ctypeRNum  t323 = t132*t165*t247*(1.1E+1/2.4E+1);
    ctypeRNum  t324 = t133*t166*t247*(1.1E+1/2.4E+1);
    ctypeRNum  t325 = t134*t167*t249*(1.1E+1/2.4E+1);
    ctypeRNum  t326 = t135*t168*t249*(1.1E+1/2.4E+1);
    ctypeRNum  t327 = t136*t169*t249*(1.1E+1/2.4E+1);
    ctypeRNum  t328 = t137*t170*t251*(1.1E+1/2.4E+1);
    ctypeRNum  t329 = t138*t171*t251*(1.1E+1/2.4E+1);
    ctypeRNum  t330 = t139*t172*t251*(1.1E+1/2.4E+1);
    ctypeRNum  t331 = t140*t173*t253*(1.1E+1/2.4E+1);
    ctypeRNum  t332 = t141*t174*t253*(1.1E+1/2.4E+1);
    ctypeRNum  t333 = t142*t175*t253*(1.1E+1/2.4E+1);
    ctypeRNum  t334 = t143*t176*t255*(1.1E+1/2.4E+1);
    ctypeRNum  t335 = t144*t177*t255*(1.1E+1/2.4E+1);
    ctypeRNum  t336 = t145*t178*t255*(1.1E+1/2.4E+1);
    ctypeRNum  t337 = t146*t179*t257*(1.1E+1/2.4E+1);
    ctypeRNum  t338 = t147*t180*t257*(1.1E+1/2.4E+1);
    ctypeRNum  t339 = t148*t181*t257*(1.1E+1/2.4E+1);
    ctypeRNum  t340 = t223+t280;
    ctypeRNum  t341 = t224+t281;
    ctypeRNum  t342 = t225+t282;
    ctypeRNum  t343 = t270+t310+2.0E+1;
    ctypeRNum  t344 = t270+t311+2.0E+1;
    ctypeRNum  t345 = t270+t312+2.0E+1;
    ctypeRNum  t346 = t271+t313+2.0E+1;
    ctypeRNum  t347 = t271+t314+2.0E+1;
    ctypeRNum  t348 = t271+t315+2.0E+1;
    ctypeRNum  t349 = t272+t316+2.0E+1;
    ctypeRNum  t350 = t272+t317+2.0E+1;
    ctypeRNum  t351 = t272+t318+2.0E+1;
    ctypeRNum  t352 = t273+t319+2.0E+1;
    ctypeRNum  t353 = t273+t320+2.0E+1;
    ctypeRNum  t354 = t273+t321+2.0E+1;
    ctypeRNum  t355 = t274+t322+2.0E+1;
    ctypeRNum  t356 = t274+t323+2.0E+1;
    ctypeRNum  t357 = t274+t324+2.0E+1;
    ctypeRNum  t358 = t275+t325+2.0E+1;
    ctypeRNum  t359 = t275+t326+2.0E+1;
    ctypeRNum  t360 = t275+t327+2.0E+1;
    ctypeRNum  t361 = t276+t328+2.0E+1;
    ctypeRNum  t362 = t276+t329+2.0E+1;
    ctypeRNum  t363 = t276+t330+2.0E+1;
    ctypeRNum  t364 = t277+t331+2.0E+1;
    ctypeRNum  t365 = t277+t332+2.0E+1;
    ctypeRNum  t366 = t277+t333+2.0E+1;
    ctypeRNum  t367 = t278+t334+2.0E+1;
    ctypeRNum  t368 = t278+t335+2.0E+1;
    ctypeRNum  t369 = t278+t336+2.0E+1;
    ctypeRNum  t370 = t279+t337+2.0E+1;
    ctypeRNum  t371 = t279+t338+2.0E+1;
    ctypeRNum  t372 = t279+t339+2.0E+1;
    out[0] = t280*vec[10]+t281*vec[11]-t340*vec[4]-t341*vec[5]+t343*vec[9]-vec[3]*(t222+t343+t38*t220*(1.1E+1/1.2E+1)+2.0E+1);
    out[1] = t280*vec[9]+t282*vec[11]-t340*vec[3]-t342*vec[5]+t344*vec[10]-vec[4]*(t222+t344+t39*t220*(1.1E+1/1.2E+1)+2.0E+1);
    out[2] = t281*vec[9]+t282*vec[10]-t341*vec[3]-t342*vec[4]+t345*vec[11]-vec[5]*(t222+t345+t40*t220*(1.1E+1/1.2E+1)+2.0E+1);
    out[3] = vec[0];
    out[4] = vec[1];
    out[5] = vec[2];
    out[6] = vec[9]*-4.0E+1+t260*vec[9]+t261*vec[9]+t280*vec[4]+t281*vec[5]+t283*vec[16]+t284*vec[17]+t343*vec[3]+t346*vec[15]-vec[10]*(t280+t123*t155*t241*(1.1E+1/2.4E+1))-vec[11]*(t281+t124*t155*t241*(1.1E+1/2.4E+1))-t119*t152*t239*vec[9]*(1.1E+1/2.4E+1)-t122*t155*t241*vec[9]*(1.1E+1/2.4E+1);
    out[7] = vec[10]*-4.0E+1+t260*vec[10]+t261*vec[10]+t280*vec[3]+t282*vec[5]+t283*vec[15]+t285*vec[17]+t344*vec[4]+t347*vec[16]-vec[9]*(t280+t122*t156*t241*(1.1E+1/2.4E+1))-vec[11]*(t282+t124*t156*t241*(1.1E+1/2.4E+1))-t120*t153*t239*vec[10]*(1.1E+1/2.4E+1)-t123*t156*t241*vec[10]*(1.1E+1/2.4E+1);
    out[8] = vec[11]*-4.0E+1+t260*vec[11]+t261*vec[11]+t281*vec[3]+t282*vec[4]+t284*vec[15]+t285*vec[16]+t345*vec[5]+t348*vec[17]-vec[9]*(t281+t122*t157*t241*(1.1E+1/2.4E+1))-vec[10]*(t282+t123*t157*t241*(1.1E+1/2.4E+1))-t121*t154*t239*vec[11]*(1.1E+1/2.4E+1)-t124*t157*t241*vec[11]*(1.1E+1/2.4E+1);
    out[9] = vec[6];
    out[10] = vec[7];
    out[11] = vec[8];
    out[12] = vec[15]*-4.0E+1+t261*vec[15]+t262*vec[15]+t283*vec[10]+t284*vec[11]+t286*vec[22]+t287*vec[23]+t346*vec[9]+t349*vec[21]-vec[16]*(t283+t126*t158*t243*(1.1E+1/2.4E+1))-vec[17]*(t284+t127*t158*t243*(1.1E+1/2.4E+1))-t122*t155*t241*vec[15]*(1.1E+1/2.4E+1)-t125*t158*t243*vec[15]*(1.1E+1/2.4E+1);
    out[13] = vec[16]*-4.0E+1+t261*vec[16]+t262*vec[16]+t283*vec[9]+t285*vec[11]+t286*vec[21]+t288*vec[23]+t347*vec[10]+t350*vec[22]-vec[15]*(t283+t125*t159*t243*(1.1E+1/2.4E+1))-vec[17]*(t285+t127*t159*t243*(1.1E+1/2.4E+1))-t123*t156*t241*vec[16]*(1.1E+1/2.4E+1)-t126*t159*t243*vec[16]*(1.1E+1/2.4E+1);
    out[14] = vec[17]*-4.0E+1+t261*vec[17]+t262*vec[17]+t284*vec[9]+t285*vec[10]+t287*vec[21]+t288*vec[22]+t348*vec[11]+t351*vec[23]-vec[15]*(t284+t125*t160*t243*(1.1E+1/2.4E+1))-vec[16]*(t285+t126*t160*t243*(1.1E+1/2.4E+1))-t124*t157*t241*vec[17]*(1.1E+1/2.4E+1)-t127*t160*t243*vec[17]*(1.1E+1/2.4E+1);
    out[15] = vec[12];
    out[16] = vec[13];
    out[17] = vec[14];
    out[18] = vec[21]*-4.0E+1+t262*vec[21]+t263*vec[21]+t286*vec[16]+t287*vec[17]+t289*vec[28]+t290*vec[29]+t349*vec[15]+t352*vec[27]-vec[22]*(t286+t129*t161*t245*(1.1E+1/2.4E+1))-vec[23]*(t287+t130*t161*t245*(1.1E+1/2.4E+1))-t125*t158*t243*vec[21]*(1.1E+1/2.4E+1)-t128*t161*t245*vec[21]*(1.1E+1/2.4E+1);
    out[19] = vec[22]*-4.0E+1+t262*vec[22]+t263*vec[22]+t286*vec[15]+t288*vec[17]+t289*vec[27]+t291*vec[29]+t350*vec[16]+t353*vec[28]-vec[21]*(t286+t128*t162*t245*(1.1E+1/2.4E+1))-vec[23]*(t288+t130*t162*t245*(1.1E+1/2.4E+1))-t126*t159*t243*vec[22]*(1.1E+1/2.4E+1)-t129*t162*t245*vec[22]*(1.1E+1/2.4E+1);
    out[20] = vec[23]*-4.0E+1+t262*vec[23]+t263*vec[23]+t287*vec[15]+t288*vec[16]+t290*vec[27]+t291*vec[28]+t351*vec[17]+t354*vec[29]-vec[21]*(t287+t128*t163*t245*(1.1E+1/2.4E+1))-vec[22]*(t288+t129*t163*t245*(1.1E+1/2.4E+1))-t127*t160*t243*vec[23]*(1.1E+1/2.4E+1)-t130*t163*t245*vec[23]*(1.1E+1/2.4E+1);
    out[21] = vec[18];
    out[22] = vec[19];
    out[23] = vec[20];
    out[24] = vec[27]*-4.0E+1+t263*vec[27]+t264*vec[27]+t289*vec[22]+t290*vec[23]+t292*vec[34]+t293*vec[35]+t352*vec[21]+t355*vec[33]-vec[28]*(t289+t132*t164*t247*(1.1E+1/2.4E+1))-vec[29]*(t290+t133*t164*t247*(1.1E+1/2.4E+1))-t128*t161*t245*vec[27]*(1.1E+1/2.4E+1)-t131*t164*t247*vec[27]*(1.1E+1/2.4E+1);
    out[25] = vec[28]*-4.0E+1+t263*vec[28]+t264*vec[28]+t289*vec[21]+t291*vec[23]+t292*vec[33]+t294*vec[35]+t353*vec[22]+t356*vec[34]-vec[27]*(t289+t131*t165*t247*(1.1E+1/2.4E+1))-vec[29]*(t291+t133*t165*t247*(1.1E+1/2.4E+1))-t129*t162*t245*vec[28]*(1.1E+1/2.4E+1)-t132*t165*t247*vec[28]*(1.1E+1/2.4E+1);
    out[26] = vec[29]*-4.0E+1+t263*vec[29]+t264*vec[29]+t290*vec[21]+t291*vec[22]+t293*vec[33]+t294*vec[34]+t354*vec[23]+t357*vec[35]-vec[27]*(t290+t131*t166*t247*(1.1E+1/2.4E+1))-vec[28]*(t291+t132*t166*t247*(1.1E+1/2.4E+1))-t130*t163*t245*vec[29]*(1.1E+1/2.4E+1)-t133*t166*t247*vec[29]*(1.1E+1/2.4E+1);
    out[27] = vec[24];
    out[28] = vec[25];
    out[29] = vec[26];
    out[30] = vec[33]*-4.0E+1+t264*vec[33]+t265*vec[33]+t292*vec[28]+t293*vec[29]+t295*vec[40]+t296*vec[41]+t355*vec[27]+t358*vec[39]-vec[34]*(t292+t135*t167*t249*(1.1E+1/2.4E+1))-vec[35]*(t293+t136*t167*t249*(1.1E+1/2.4E+1))-t131*t164*t247*vec[33]*(1.1E+1/2.4E+1)-t134*t167*t249*vec[33]*(1.1E+1/2.4E+1);
    out[31] = vec[34]*-4.0E+1+t264*vec[34]+t265*vec[34]+t292*vec[27]+t294*vec[29]+t295*vec[39]+t297*vec[41]+t356*vec[28]+t359*vec[40]-vec[33]*(t292+t134*t168*t249*(1.1E+1/2.4E+1))-vec[35]*(t294+t136*t168*t249*(1.1E+1/2.4E+1))-t132*t165*t247*vec[34]*(1.1E+1/2.4E+1)-t135*t168*t249*vec[34]*(1.1E+1/2.4E+1);
    out[32] = vec[35]*-4.0E+1+t264*vec[35]+t265*vec[35]+t293*vec[27]+t294*vec[28]+t296*vec[39]+t297*vec[40]+t357*vec[29]+t360*vec[41]-vec[33]*(t293+t134*t169*t249*(1.1E+1/2.4E+1))-vec[34]*(t294+t135*t169*t249*(1.1E+1/2.4E+1))-t133*t166*t247*vec[35]*(1.1E+1/2.4E+1)-t136*t169*t249*vec[35]*(1.1E+1/2.4E+1);
    out[33] = vec[30];
    out[34] = vec[31];
    out[35] = vec[32];
    out[36] = vec[39]*-4.0E+1+t265*vec[39]+t266*vec[39]+t295*vec[34]+t296*vec[35]+t298*vec[46]+t299*vec[47]+t358*vec[33]+t361*vec[45]-vec[40]*(t295+t138*t170*t251*(1.1E+1/2.4E+1))-vec[41]*(t296+t139*t170*t251*(1.1E+1/2.4E+1))-t134*t167*t249*vec[39]*(1.1E+1/2.4E+1)-t137*t170*t251*vec[39]*(1.1E+1/2.4E+1);
    out[37] = vec[40]*-4.0E+1+t265*vec[40]+t266*vec[40]+t295*vec[33]+t297*vec[35]+t298*vec[45]+t300*vec[47]+t359*vec[34]+t362*vec[46]-vec[39]*(t295+t137*t171*t251*(1.1E+1/2.4E+1))-vec[41]*(t297+t139*t171*t251*(1.1E+1/2.4E+1))-t135*t168*t249*vec[40]*(1.1E+1/2.4E+1)-t138*t171*t251*vec[40]*(1.1E+1/2.4E+1);
    out[38] = vec[41]*-4.0E+1+t265*vec[41]+t266*vec[41]+t296*vec[33]+t297*vec[34]+t299*vec[45]+t300*vec[46]+t360*vec[35]+t363*vec[47]-vec[39]*(t296+t137*t172*t251*(1.1E+1/2.4E+1))-vec[40]*(t297+t138*t172*t251*(1.1E+1/2.4E+1))-t136*t169*t249*vec[41]*(1.1E+1/2.4E+1)-t139*t172*t251*vec[41]*(1.1E+1/2.4E+1);
    out[39] = vec[36];
    out[40] = vec[37];
    out[41] = vec[38];
    out[42] = vec[45]*-4.0E+1+t266*vec[45]+t267*vec[45]+t298*vec[40]+t299*vec[41]+t301*vec[52]+t302*vec[53]+t361*vec[39]+t364*vec[51]-vec[46]*(t298+t141*t173*t253*(1.1E+1/2.4E+1))-vec[47]*(t299+t142*t173*t253*(1.1E+1/2.4E+1))-t137*t170*t251*vec[45]*(1.1E+1/2.4E+1)-t140*t173*t253*vec[45]*(1.1E+1/2.4E+1);
    out[43] = vec[46]*-4.0E+1+t266*vec[46]+t267*vec[46]+t298*vec[39]+t300*vec[41]+t301*vec[51]+t303*vec[53]+t362*vec[40]+t365*vec[52]-vec[45]*(t298+t140*t174*t253*(1.1E+1/2.4E+1))-vec[47]*(t300+t142*t174*t253*(1.1E+1/2.4E+1))-t138*t171*t251*vec[46]*(1.1E+1/2.4E+1)-t141*t174*t253*vec[46]*(1.1E+1/2.4E+1);
    out[44] = vec[47]*-4.0E+1+t266*vec[47]+t267*vec[47]+t299*vec[39]+t300*vec[40]+t302*vec[51]+t303*vec[52]+t363*vec[41]+t366*vec[53]-vec[45]*(t299+t140*t175*t253*(1.1E+1/2.4E+1))-vec[46]*(t300+t141*t175*t253*(1.1E+1/2.4E+1))-t139*t172*t251*vec[47]*(1.1E+1/2.4E+1)-t142*t175*t253*vec[47]*(1.1E+1/2.4E+1);
    out[45] = vec[42];
    out[46] = vec[43];
    out[47] = vec[44];
    out[48] = vec[51]*-4.0E+1+t267*vec[51]+t268*vec[51]+t301*vec[46]+t302*vec[47]+t304*vec[58]+t305*vec[59]+t364*vec[45]+t367*vec[57]-vec[52]*(t301+t144*t176*t255*(1.1E+1/2.4E+1))-vec[53]*(t302+t145*t176*t255*(1.1E+1/2.4E+1))-t140*t173*t253*vec[51]*(1.1E+1/2.4E+1)-t143*t176*t255*vec[51]*(1.1E+1/2.4E+1);
    out[49] = vec[52]*-4.0E+1+t267*vec[52]+t268*vec[52]+t301*vec[45]+t303*vec[47]+t304*vec[57]+t306*vec[59]+t365*vec[46]+t368*vec[58]-vec[51]*(t301+t143*t177*t255*(1.1E+1/2.4E+1))-vec[53]*(t303+t145*t177*t255*(1.1E+1/2.4E+1))-t141*t174*t253*vec[52]*(1.1E+1/2.4E+1)-t144*t177*t255*vec[52]*(1.1E+1/2.4E+1);
    out[50] = vec[53]*-4.0E+1+t267*vec[53]+t268*vec[53]+t302*vec[45]+t303*vec[46]+t305*vec[57]+t306*vec[58]+t366*vec[47]+t369*vec[59]-vec[51]*(t302+t143*t178*t255*(1.1E+1/2.4E+1))-vec[52]*(t303+t144*t178*t255*(1.1E+1/2.4E+1))-t142*t175*t253*vec[53]*(1.1E+1/2.4E+1)-t145*t178*t255*vec[53]*(1.1E+1/2.4E+1);
    out[51] = vec[48];
    out[52] = vec[49];
    out[53] = vec[50];
    out[54] = vec[57]*-4.0E+1+t268*vec[57]+t269*vec[57]+t304*vec[52]+t305*vec[53]+t307*vec[64]+t308*vec[65]+t367*vec[51]+t370*vec[63]-vec[58]*(t304+t147*t179*t257*(1.1E+1/2.4E+1))-vec[59]*(t305+t148*t179*t257*(1.1E+1/2.4E+1))-t143*t176*t255*vec[57]*(1.1E+1/2.4E+1)-t146*t179*t257*vec[57]*(1.1E+1/2.4E+1);
    out[55] = vec[58]*-4.0E+1+t268*vec[58]+t269*vec[58]+t304*vec[51]+t306*vec[53]+t307*vec[63]+t309*vec[65]+t368*vec[52]+t371*vec[64]-vec[57]*(t304+t146*t180*t257*(1.1E+1/2.4E+1))-vec[59]*(t306+t148*t180*t257*(1.1E+1/2.4E+1))-t144*t177*t255*vec[58]*(1.1E+1/2.4E+1)-t147*t180*t257*vec[58]*(1.1E+1/2.4E+1);
    out[56] = vec[59]*-4.0E+1+t268*vec[59]+t269*vec[59]+t305*vec[51]+t306*vec[52]+t308*vec[63]+t309*vec[64]+t369*vec[53]+t372*vec[65]-vec[57]*(t305+t146*t181*t257*(1.1E+1/2.4E+1))-vec[58]*(t306+t147*t181*t257*(1.1E+1/2.4E+1))-t145*t178*t255*vec[59]*(1.1E+1/2.4E+1)-t148*t181*t257*vec[59]*(1.1E+1/2.4E+1);
    out[57] = vec[54];
    out[58] = vec[55];
    out[59] = vec[56];
    out[60] = vec[63]*-4.0E+1+t258*vec[63]*(1.1E+1/1.2E+1)+t269*vec[63]+t307*vec[58]+t308*vec[59]+t370*vec[57]-vec[64]*(t307+t150*t182*t259*(1.1E+1/2.4E+1))-vec[65]*(t308+t151*t182*t259*(1.1E+1/2.4E+1))-t146*t179*t257*vec[63]*(1.1E+1/2.4E+1)-t149*t182*t259*vec[63]*(1.1E+1/2.4E+1);
    out[61] = vec[64]*-4.0E+1+t258*vec[64]*(1.1E+1/1.2E+1)+t269*vec[64]+t307*vec[57]+t309*vec[59]+t371*vec[58]-vec[63]*(t307+t149*t183*t259*(1.1E+1/2.4E+1))-vec[65]*(t309+t151*t183*t259*(1.1E+1/2.4E+1))-t147*t180*t257*vec[64]*(1.1E+1/2.4E+1)-t150*t183*t259*vec[64]*(1.1E+1/2.4E+1);
    out[62] = vec[65]*-4.0E+1+t258*vec[65]*(1.1E+1/1.2E+1)+t269*vec[65]+t308*vec[57]+t309*vec[58]+t372*vec[59]-vec[63]*(t308+t149*t184*t259*(1.1E+1/2.4E+1))-vec[64]*(t309+t150*t184*t259*(1.1E+1/2.4E+1))-t148*t181*t257*vec[65]*(1.1E+1/2.4E+1)-t151*t184*t259*vec[65]*(1.1E+1/2.4E+1);
    out[63] = vec[60];
    out[64] = vec[61];
    out[65] = vec[62];
    out[66] = (t259*(t41*vec[63]*1.1E+1+t237*vec[63]*2.4E+2+t111*x[60]+t112*x[60]-vec[63]*(t41+t216+t217)*1.1E+1-vec[64]*x[60]*x[67]*1.1E+1-vec[65]*x[60]*x[68]*1.1E+1))/1.2E+1-(t259*x[66]*(t111+t112+t117+t118))/1.2E+1;
    out[67] = (t259*(t42*vec[64]*1.1E+1+t237*vec[64]*2.4E+2+t110*x[61]+t112*x[61]-vec[64]*(t42+t215+t217)*1.1E+1-vec[63]*x[61]*x[66]*1.1E+1-vec[65]*x[61]*x[68]*1.1E+1))/1.2E+1-(t259*x[67]*(t110+t112+t116+t118))/1.2E+1;
    out[68] = (t259*(t43*vec[65]*1.1E+1+t237*vec[65]*2.4E+2+t110*x[62]+t111*x[62]-vec[65]*(t43+t215+t216)*1.1E+1-vec[63]*x[62]*x[66]*1.1E+1-vec[64]*x[62]*x[67]*1.1E+1))/1.2E+1-(t259*x[68]*(t110+t111+t116+t117))/1.2E+1;
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Chain_12_ProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
    out[0] = vec[66];
    out[1] = vec[67];
    out[2] = vec[68];
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void Chain_12_ProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void Chain_12_ProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    

    out[0] = pCost_[2]*(u[0]*u[0])+pCost_[2]*(u[1]*u[1])+pCost_[2]*(u[2]*u[2])+pCost_[1]*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+pCost_[1]*(x[9]*x[9]+x[10]*x[10]+x[11]*x[11])+pCost_[1]*(x[15]*x[15]+x[16]*x[16]+x[17]*x[17])+pCost_[1]*(x[21]*x[21]+x[22]*x[22]+x[23]*x[23])+pCost_[1]*(x[27]*x[27]+x[28]*x[28]+x[29]*x[29])+pCost_[1]*(x[33]*x[33]+x[34]*x[34]+x[35]*x[35])+pCost_[1]*(x[39]*x[39]+x[40]*x[40]+x[41]*x[41])+pCost_[1]*(x[45]*x[45]+x[46]*x[46]+x[47]*x[47])+pCost_[1]*(x[51]*x[51]+x[52]*x[52]+x[53]*x[53])+pCost_[1]*(x[57]*x[57]+x[58]*x[58]+x[59]*x[59])+pCost_[1]*(x[63]*x[63]+x[64]*x[64]+x[65]*x[65])+pCost_[0]*(POW2(x[66]-6.0)+x[67]*x[67]+x[68]*x[68]);
}


/** Gradient dl/dx **/
void Chain_12_ProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    

    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = pCost_[1]*x[3]*2.0;
    out[4] = pCost_[1]*x[4]*2.0;
    out[5] = pCost_[1]*x[5]*2.0;
    out[6] = 0.0;
    out[7] = 0.0;
    out[8] = 0.0;
    out[9] = pCost_[1]*x[9]*2.0;
    out[10] = pCost_[1]*x[10]*2.0;
    out[11] = pCost_[1]*x[11]*2.0;
    out[12] = 0.0;
    out[13] = 0.0;
    out[14] = 0.0;
    out[15] = pCost_[1]*x[15]*2.0;
    out[16] = pCost_[1]*x[16]*2.0;
    out[17] = pCost_[1]*x[17]*2.0;
    out[18] = 0.0;
    out[19] = 0.0;
    out[20] = 0.0;
    out[21] = pCost_[1]*x[21]*2.0;
    out[22] = pCost_[1]*x[22]*2.0;
    out[23] = pCost_[1]*x[23]*2.0;
    out[24] = 0.0;
    out[25] = 0.0;
    out[26] = 0.0;
    out[27] = pCost_[1]*x[27]*2.0;
    out[28] = pCost_[1]*x[28]*2.0;
    out[29] = pCost_[1]*x[29]*2.0;
    out[30] = 0.0;
    out[31] = 0.0;
    out[32] = 0.0;
    out[33] = pCost_[1]*x[33]*2.0;
    out[34] = pCost_[1]*x[34]*2.0;
    out[35] = pCost_[1]*x[35]*2.0;
    out[36] = 0.0;
    out[37] = 0.0;
    out[38] = 0.0;
    out[39] = pCost_[1]*x[39]*2.0;
    out[40] = pCost_[1]*x[40]*2.0;
    out[41] = pCost_[1]*x[41]*2.0;
    out[42] = 0.0;
    out[43] = 0.0;
    out[44] = 0.0;
    out[45] = pCost_[1]*x[45]*2.0;
    out[46] = pCost_[1]*x[46]*2.0;
    out[47] = pCost_[1]*x[47]*2.0;
    out[48] = 0.0;
    out[49] = 0.0;
    out[50] = 0.0;
    out[51] = pCost_[1]*x[51]*2.0;
    out[52] = pCost_[1]*x[52]*2.0;
    out[53] = pCost_[1]*x[53]*2.0;
    out[54] = 0.0;
    out[55] = 0.0;
    out[56] = 0.0;
    out[57] = pCost_[1]*x[57]*2.0;
    out[58] = pCost_[1]*x[58]*2.0;
    out[59] = pCost_[1]*x[59]*2.0;
    out[60] = 0.0;
    out[61] = 0.0;
    out[62] = 0.0;
    out[63] = pCost_[1]*x[63]*2.0;
    out[64] = pCost_[1]*x[64]*2.0;
    out[65] = pCost_[1]*x[65]*2.0;
    out[66] = pCost_[0]*(x[66]-6.0)*2.0;
    out[67] = pCost_[0]*x[67]*2.0;
    out[68] = pCost_[0]*x[68]*2.0;
}


/** Gradient dl/du **/
void Chain_12_ProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    

    out[0] = pCost_[2]*u[0]*2.0;
    out[1] = pCost_[2]*u[1]*2.0;
    out[2] = pCost_[2]*u[2]*2.0;
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Chain_12_ProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 0.0;
}


/** Gradient dV/dx **/
void Chain_12_ProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
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
}


/** Gradient dV/dT **/
void Chain_12_ProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void Chain_12_ProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Chain_12_ProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Chain_12_ProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void Chain_12_ProblemDescription::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p)
{
}


/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void Chain_12_ProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
{
}


/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void Chain_12_ProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
{
}