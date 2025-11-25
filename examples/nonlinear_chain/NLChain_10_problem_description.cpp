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
* This probfct-file describes the nonlinear chain problem with 10 chain elements from
* Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
* Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
* Process Control 22(3), 540-550 (2012)
*
*/

#include "NLChain_10_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

using namespace grampc;


Chain_10_ProblemDescription::Chain_10_ProblemDescription(const std::vector<typeRNum>& pCost)
: pCost_(pCost)
{
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void Chain_10_ProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 57;
    *Nu  = 3;
    *Np  = 0;
    *Nh  = 0;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void Chain_10_ProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
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
    ctypeRNum  t32 = t5+x[0];
    ctypeRNum  t33 = t6+x[1];
    ctypeRNum  t34 = t7+x[2];
    ctypeRNum  t35 = t8+x[6];
    ctypeRNum  t36 = t9+x[7];
    ctypeRNum  t37 = t10+x[8];
    ctypeRNum  t38 = t11+x[12];
    ctypeRNum  t39 = t12+x[13];
    ctypeRNum  t40 = t13+x[14];
    ctypeRNum  t41 = t14+x[18];
    ctypeRNum  t42 = t15+x[19];
    ctypeRNum  t43 = t16+x[20];
    ctypeRNum  t44 = t17+x[24];
    ctypeRNum  t45 = t18+x[25];
    ctypeRNum  t46 = t19+x[26];
    ctypeRNum  t47 = t20+x[30];
    ctypeRNum  t48 = t21+x[31];
    ctypeRNum  t49 = t22+x[32];
    ctypeRNum  t50 = t23+x[36];
    ctypeRNum  t51 = t24+x[37];
    ctypeRNum  t52 = t25+x[38];
    ctypeRNum  t53 = t26+x[42];
    ctypeRNum  t54 = t27+x[43];
    ctypeRNum  t55 = t28+x[44];
    ctypeRNum  t56 = t29+x[48];
    ctypeRNum  t57 = t30+x[49];
    ctypeRNum  t58 = t31+x[50];
    ctypeRNum  t86 = t2+t3+t4;
    ctypeRNum  t59 = t32*t32;
    ctypeRNum  t60 = t33*t33;
    ctypeRNum  t61 = t34*t34;
    ctypeRNum  t62 = t35*t35;
    ctypeRNum  t63 = t36*t36;
    ctypeRNum  t64 = t37*t37;
    ctypeRNum  t65 = t38*t38;
    ctypeRNum  t66 = t39*t39;
    ctypeRNum  t67 = t40*t40;
    ctypeRNum  t68 = t41*t41;
    ctypeRNum  t69 = t42*t42;
    ctypeRNum  t70 = t43*t43;
    ctypeRNum  t71 = t44*t44;
    ctypeRNum  t72 = t45*t45;
    ctypeRNum  t73 = t46*t46;
    ctypeRNum  t74 = t47*t47;
    ctypeRNum  t75 = t48*t48;
    ctypeRNum  t76 = t49*t49;
    ctypeRNum  t77 = t50*t50;
    ctypeRNum  t78 = t51*t51;
    ctypeRNum  t79 = t52*t52;
    ctypeRNum  t80 = t53*t53;
    ctypeRNum  t81 = t54*t54;
    ctypeRNum  t82 = t55*t55;
    ctypeRNum  t83 = t56*t56;
    ctypeRNum  t84 = t57*t57;
    ctypeRNum  t85 = t58*t58;
    ctypeRNum  t87 = 1.0/sqrt(t86);
    ctypeRNum  t88 = t87*(1.1E+1/4.0E+2);
    ctypeRNum  t89 = t59+t60+t61;
    ctypeRNum  t90 = t62+t63+t64;
    ctypeRNum  t91 = t65+t66+t67;
    ctypeRNum  t92 = t68+t69+t70;
    ctypeRNum  t93 = t71+t72+t73;
    ctypeRNum  t94 = t74+t75+t76;
    ctypeRNum  t95 = t77+t78+t79;
    ctypeRNum  t96 = t80+t81+t82;
    ctypeRNum  t97 = t83+t84+t85;
    ctypeRNum  t98 = t88-1.0/2.0;
    ctypeRNum  t99 = 1.0/sqrt(t89);
    ctypeRNum  t100 = 1.0/sqrt(t90);
    ctypeRNum  t101 = 1.0/sqrt(t91);
    ctypeRNum  t102 = 1.0/sqrt(t92);
    ctypeRNum  t103 = 1.0/sqrt(t93);
    ctypeRNum  t104 = 1.0/sqrt(t94);
    ctypeRNum  t105 = 1.0/sqrt(t95);
    ctypeRNum  t106 = 1.0/sqrt(t96);
    ctypeRNum  t107 = 1.0/sqrt(t97);
    ctypeRNum  t108 = t99*(1.1E+1/4.0E+2);
    ctypeRNum  t109 = t100*(1.1E+1/4.0E+2);
    ctypeRNum  t110 = t101*(1.1E+1/4.0E+2);
    ctypeRNum  t111 = t102*(1.1E+1/4.0E+2);
    ctypeRNum  t112 = t103*(1.1E+1/4.0E+2);
    ctypeRNum  t113 = t104*(1.1E+1/4.0E+2);
    ctypeRNum  t114 = t105*(1.1E+1/4.0E+2);
    ctypeRNum  t115 = t106*(1.1E+1/4.0E+2);
    ctypeRNum  t116 = t107*(1.1E+1/4.0E+2);
    ctypeRNum  t117 = t108-1.0/2.0;
    ctypeRNum  t118 = t109-1.0/2.0;
    ctypeRNum  t119 = t110-1.0/2.0;
    ctypeRNum  t120 = t111-1.0/2.0;
    ctypeRNum  t121 = t112-1.0/2.0;
    ctypeRNum  t122 = t113-1.0/2.0;
    ctypeRNum  t123 = t114-1.0/2.0;
    ctypeRNum  t124 = t115-1.0/2.0;
    ctypeRNum  t125 = t116-1.0/2.0;
    ctypeRNum  t126 = t32*t117*(1.0E+2/3.0);
    ctypeRNum  t127 = t33*t117*(1.0E+2/3.0);
    ctypeRNum  t128 = t34*t117*(1.0E+2/3.0);
    ctypeRNum  t129 = t35*t118*(1.0E+2/3.0);
    ctypeRNum  t130 = t36*t118*(1.0E+2/3.0);
    ctypeRNum  t131 = t37*t118*(1.0E+2/3.0);
    ctypeRNum  t132 = t38*t119*(1.0E+2/3.0);
    ctypeRNum  t133 = t39*t119*(1.0E+2/3.0);
    ctypeRNum  t134 = t40*t119*(1.0E+2/3.0);
    ctypeRNum  t135 = t41*t120*(1.0E+2/3.0);
    ctypeRNum  t136 = t42*t120*(1.0E+2/3.0);
    ctypeRNum  t137 = t43*t120*(1.0E+2/3.0);
    ctypeRNum  t138 = t44*t121*(1.0E+2/3.0);
    ctypeRNum  t139 = t45*t121*(1.0E+2/3.0);
    ctypeRNum  t140 = t46*t121*(1.0E+2/3.0);
    ctypeRNum  t141 = t47*t122*(1.0E+2/3.0);
    ctypeRNum  t142 = t48*t122*(1.0E+2/3.0);
    ctypeRNum  t143 = t49*t122*(1.0E+2/3.0);
    ctypeRNum  t144 = t50*t123*(1.0E+2/3.0);
    ctypeRNum  t145 = t51*t123*(1.0E+2/3.0);
    ctypeRNum  t146 = t52*t123*(1.0E+2/3.0);
    ctypeRNum  t147 = t53*t124*(1.0E+2/3.0);
    ctypeRNum  t148 = t54*t124*(1.0E+2/3.0);
    ctypeRNum  t149 = t55*t124*(1.0E+2/3.0);
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = t126+t98*x[0]*(1.0E+2/3.0);
    out[4] = t127+t98*x[1]*(1.0E+2/3.0);
    out[5] = t128+t98*x[2]*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[6] = x[9];
    out[7] = x[10];
    out[8] = x[11];
    out[9] = -t126+t129;
    out[10] = -t127+t130;
    out[11] = -t128+t131-9.81E+2/1.0E+2;
    out[12] = x[15];
    out[13] = x[16];
    out[14] = x[17];
    out[15] = -t129+t132;
    out[16] = -t130+t133;
    out[17] = -t131+t134-9.81E+2/1.0E+2;
    out[18] = x[21];
    out[19] = x[22];
    out[20] = x[23];
    out[21] = -t132+t135;
    out[22] = -t133+t136;
    out[23] = -t134+t137-9.81E+2/1.0E+2;
    out[24] = x[27];
    out[25] = x[28];
    out[26] = x[29];
    out[27] = -t135+t138;
    out[28] = -t136+t139;
    out[29] = -t137+t140-9.81E+2/1.0E+2;
    out[30] = x[33];
    out[31] = x[34];
    out[32] = x[35];
    out[33] = -t138+t141;
    out[34] = -t139+t142;
    out[35] = -t140+t143-9.81E+2/1.0E+2;
    out[36] = x[39];
    out[37] = x[40];
    out[38] = x[41];
    out[39] = -t141+t144;
    out[40] = -t142+t145;
    out[41] = -t143+t146-9.81E+2/1.0E+2;
    out[42] = x[45];
    out[43] = x[46];
    out[44] = x[47];
    out[45] = -t144+t147;
    out[46] = -t145+t148;
    out[47] = -t146+t149-9.81E+2/1.0E+2;
    out[48] = x[51];
    out[49] = x[52];
    out[50] = x[53];
    out[51] = -t147+t56*t125*(1.0E+2/3.0);
    out[52] = -t148+t57*t125*(1.0E+2/3.0);
    out[53] = -t149+t58*t125*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[54] = u[0];
    out[55] = u[1];
    out[56] = u[2];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Chain_10_ProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
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
    ctypeRNum  t32 = x[0]*x[0];
    ctypeRNum  t33 = x[1]*x[1];
    ctypeRNum  t34 = x[2]*x[2];
    ctypeRNum  t35 = x[48]*x[48];
    ctypeRNum  t36 = x[49]*x[49];
    ctypeRNum  t37 = x[50]*x[50];
    ctypeRNum  t38 = -x[6];
    ctypeRNum  t40 = -x[7];
    ctypeRNum  t42 = -x[8];
    ctypeRNum  t44 = -x[12];
    ctypeRNum  t46 = -x[13];
    ctypeRNum  t48 = -x[14];
    ctypeRNum  t50 = -x[18];
    ctypeRNum  t52 = -x[19];
    ctypeRNum  t54 = -x[20];
    ctypeRNum  t56 = -x[24];
    ctypeRNum  t58 = -x[25];
    ctypeRNum  t60 = -x[26];
    ctypeRNum  t62 = -x[30];
    ctypeRNum  t64 = -x[31];
    ctypeRNum  t66 = -x[32];
    ctypeRNum  t68 = -x[36];
    ctypeRNum  t70 = -x[37];
    ctypeRNum  t72 = -x[38];
    ctypeRNum  t74 = -x[42];
    ctypeRNum  t76 = -x[43];
    ctypeRNum  t78 = -x[44];
    ctypeRNum  t80 = -x[48];
    ctypeRNum  t82 = -x[49];
    ctypeRNum  t84 = -x[50];
    ctypeRNum  t86 = -x[54];
    ctypeRNum  t88 = -x[55];
    ctypeRNum  t90 = -x[56];
    ctypeRNum  t92 = vec[51]*x[48]*1.1E+1;
    ctypeRNum  t93 = vec[52]*x[49]*1.1E+1;
    ctypeRNum  t94 = vec[53]*x[50]*1.1E+1;
    ctypeRNum  t95 = vec[51]*x[54]*1.1E+1;
    ctypeRNum  t96 = vec[52]*x[55]*1.1E+1;
    ctypeRNum  t97 = vec[53]*x[56]*1.1E+1;
    ctypeRNum  t39 = -t5;
    ctypeRNum  t41 = -t6;
    ctypeRNum  t43 = -t7;
    ctypeRNum  t45 = -t8;
    ctypeRNum  t47 = -t9;
    ctypeRNum  t49 = -t10;
    ctypeRNum  t51 = -t11;
    ctypeRNum  t53 = -t12;
    ctypeRNum  t55 = -t13;
    ctypeRNum  t57 = -t14;
    ctypeRNum  t59 = -t15;
    ctypeRNum  t61 = -t16;
    ctypeRNum  t63 = -t17;
    ctypeRNum  t65 = -t18;
    ctypeRNum  t67 = -t19;
    ctypeRNum  t69 = -t20;
    ctypeRNum  t71 = -t21;
    ctypeRNum  t73 = -t22;
    ctypeRNum  t75 = -t23;
    ctypeRNum  t77 = -t24;
    ctypeRNum  t79 = -t25;
    ctypeRNum  t81 = -t26;
    ctypeRNum  t83 = -t27;
    ctypeRNum  t85 = -t28;
    ctypeRNum  t87 = -t29;
    ctypeRNum  t89 = -t30;
    ctypeRNum  t91 = -t31;
    ctypeRNum  t98 = -t95;
    ctypeRNum  t99 = -t96;
    ctypeRNum  t100 = -t97;
    ctypeRNum  t101 = t38+x[0];
    ctypeRNum  t102 = t40+x[1];
    ctypeRNum  t103 = t42+x[2];
    ctypeRNum  t104 = t44+x[6];
    ctypeRNum  t105 = t46+x[7];
    ctypeRNum  t106 = t48+x[8];
    ctypeRNum  t107 = t50+x[12];
    ctypeRNum  t108 = t52+x[13];
    ctypeRNum  t109 = t54+x[14];
    ctypeRNum  t110 = t56+x[18];
    ctypeRNum  t111 = t58+x[19];
    ctypeRNum  t112 = t60+x[20];
    ctypeRNum  t113 = t62+x[24];
    ctypeRNum  t114 = t64+x[25];
    ctypeRNum  t115 = t66+x[26];
    ctypeRNum  t116 = t68+x[30];
    ctypeRNum  t117 = t70+x[31];
    ctypeRNum  t118 = t72+x[32];
    ctypeRNum  t119 = t74+x[36];
    ctypeRNum  t120 = t76+x[37];
    ctypeRNum  t121 = t78+x[38];
    ctypeRNum  t122 = t80+x[42];
    ctypeRNum  t123 = t82+x[43];
    ctypeRNum  t124 = t84+x[44];
    ctypeRNum  t125 = t86+x[48];
    ctypeRNum  t126 = t88+x[49];
    ctypeRNum  t127 = t90+x[50];
    ctypeRNum  t182 = t32+t33+t34;
    ctypeRNum  t128 = t2+t39;
    ctypeRNum  t129 = t3+t41;
    ctypeRNum  t130 = t4+t43;
    ctypeRNum  t131 = t5+t45;
    ctypeRNum  t132 = t6+t47;
    ctypeRNum  t133 = t7+t49;
    ctypeRNum  t134 = t8+t51;
    ctypeRNum  t135 = t9+t53;
    ctypeRNum  t136 = t10+t55;
    ctypeRNum  t137 = t11+t57;
    ctypeRNum  t138 = t12+t59;
    ctypeRNum  t139 = t13+t61;
    ctypeRNum  t140 = t14+t63;
    ctypeRNum  t141 = t15+t65;
    ctypeRNum  t142 = t16+t67;
    ctypeRNum  t143 = t17+t69;
    ctypeRNum  t144 = t18+t71;
    ctypeRNum  t145 = t19+t73;
    ctypeRNum  t146 = t20+t75;
    ctypeRNum  t147 = t21+t77;
    ctypeRNum  t148 = t22+t79;
    ctypeRNum  t149 = t23+t81;
    ctypeRNum  t150 = t24+t83;
    ctypeRNum  t151 = t25+t85;
    ctypeRNum  t152 = t26+t87;
    ctypeRNum  t153 = t27+t89;
    ctypeRNum  t154 = t28+t91;
    ctypeRNum  t155 = t101*t101;
    ctypeRNum  t156 = t102*t102;
    ctypeRNum  t157 = t103*t103;
    ctypeRNum  t158 = t104*t104;
    ctypeRNum  t159 = t105*t105;
    ctypeRNum  t160 = t106*t106;
    ctypeRNum  t161 = t107*t107;
    ctypeRNum  t162 = t108*t108;
    ctypeRNum  t163 = t109*t109;
    ctypeRNum  t164 = t110*t110;
    ctypeRNum  t165 = t111*t111;
    ctypeRNum  t166 = t112*t112;
    ctypeRNum  t167 = t113*t113;
    ctypeRNum  t168 = t114*t114;
    ctypeRNum  t169 = t115*t115;
    ctypeRNum  t170 = t116*t116;
    ctypeRNum  t171 = t117*t117;
    ctypeRNum  t172 = t118*t118;
    ctypeRNum  t173 = t119*t119;
    ctypeRNum  t174 = t120*t120;
    ctypeRNum  t175 = t121*t121;
    ctypeRNum  t176 = t122*t122;
    ctypeRNum  t177 = t123*t123;
    ctypeRNum  t178 = t124*t124;
    ctypeRNum  t179 = t125*t125;
    ctypeRNum  t180 = t126*t126;
    ctypeRNum  t181 = t127*t127;
    ctypeRNum  t183 = 1.0/sqrt(t182);
    ctypeRNum  t184 = t183*t183*t183;
    ctypeRNum  t185 = t183*(1.1E+1/1.2E+1);
    ctypeRNum  t190 = t155+t156+t157;
    ctypeRNum  t191 = t158+t159+t160;
    ctypeRNum  t192 = t161+t162+t163;
    ctypeRNum  t193 = t164+t165+t166;
    ctypeRNum  t194 = t167+t168+t169;
    ctypeRNum  t195 = t170+t171+t172;
    ctypeRNum  t196 = t173+t174+t175;
    ctypeRNum  t197 = t176+t177+t178;
    ctypeRNum  t198 = t179+t180+t181;
    ctypeRNum  t186 = -t185;
    ctypeRNum  t187 = t184*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t188 = t184*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t189 = t184*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t199 = pow(t198,3.0/2.0);
    ctypeRNum  t200 = 1.0/sqrt(t190);
    ctypeRNum  t202 = 1.0/sqrt(t191);
    ctypeRNum  t204 = 1.0/sqrt(t192);
    ctypeRNum  t206 = 1.0/sqrt(t193);
    ctypeRNum  t208 = 1.0/sqrt(t194);
    ctypeRNum  t210 = 1.0/sqrt(t195);
    ctypeRNum  t212 = 1.0/sqrt(t196);
    ctypeRNum  t214 = 1.0/sqrt(t197);
    ctypeRNum  t216 = 1.0/sqrt(t198);
    ctypeRNum  t201 = t200*t200*t200;
    ctypeRNum  t203 = t202*t202*t202;
    ctypeRNum  t205 = t204*t204*t204;
    ctypeRNum  t207 = t206*t206*t206;
    ctypeRNum  t209 = t208*t208*t208;
    ctypeRNum  t211 = t210*t210*t210;
    ctypeRNum  t213 = t212*t212*t212;
    ctypeRNum  t215 = t214*t214*t214;
    ctypeRNum  t217 = 1.0/t199;
    ctypeRNum  t218 = t200*(1.1E+1/1.2E+1);
    ctypeRNum  t219 = t202*(1.1E+1/1.2E+1);
    ctypeRNum  t220 = t204*(1.1E+1/1.2E+1);
    ctypeRNum  t221 = t206*(1.1E+1/1.2E+1);
    ctypeRNum  t222 = t208*(1.1E+1/1.2E+1);
    ctypeRNum  t223 = t210*(1.1E+1/1.2E+1);
    ctypeRNum  t224 = t212*(1.1E+1/1.2E+1);
    ctypeRNum  t225 = t214*(1.1E+1/1.2E+1);
    ctypeRNum  t226 = -t218;
    ctypeRNum  t227 = -t219;
    ctypeRNum  t228 = -t220;
    ctypeRNum  t229 = -t221;
    ctypeRNum  t230 = -t222;
    ctypeRNum  t231 = -t223;
    ctypeRNum  t232 = -t224;
    ctypeRNum  t233 = -t225;
    ctypeRNum  t234 = t101*t102*t201*(1.1E+1/1.2E+1);
    ctypeRNum  t235 = t101*t103*t201*(1.1E+1/1.2E+1);
    ctypeRNum  t236 = t102*t103*t201*(1.1E+1/1.2E+1);
    ctypeRNum  t237 = t104*t105*t203*(1.1E+1/1.2E+1);
    ctypeRNum  t238 = t104*t106*t203*(1.1E+1/1.2E+1);
    ctypeRNum  t239 = t105*t106*t203*(1.1E+1/1.2E+1);
    ctypeRNum  t240 = t107*t108*t205*(1.1E+1/1.2E+1);
    ctypeRNum  t241 = t107*t109*t205*(1.1E+1/1.2E+1);
    ctypeRNum  t242 = t108*t109*t205*(1.1E+1/1.2E+1);
    ctypeRNum  t243 = t110*t111*t207*(1.1E+1/1.2E+1);
    ctypeRNum  t244 = t110*t112*t207*(1.1E+1/1.2E+1);
    ctypeRNum  t245 = t111*t112*t207*(1.1E+1/1.2E+1);
    ctypeRNum  t246 = t113*t114*t209*(1.1E+1/1.2E+1);
    ctypeRNum  t247 = t113*t115*t209*(1.1E+1/1.2E+1);
    ctypeRNum  t248 = t114*t115*t209*(1.1E+1/1.2E+1);
    ctypeRNum  t249 = t116*t117*t211*(1.1E+1/1.2E+1);
    ctypeRNum  t250 = t116*t118*t211*(1.1E+1/1.2E+1);
    ctypeRNum  t251 = t117*t118*t211*(1.1E+1/1.2E+1);
    ctypeRNum  t252 = t119*t120*t213*(1.1E+1/1.2E+1);
    ctypeRNum  t253 = t119*t121*t213*(1.1E+1/1.2E+1);
    ctypeRNum  t254 = t120*t121*t213*(1.1E+1/1.2E+1);
    ctypeRNum  t255 = t122*t123*t215*(1.1E+1/1.2E+1);
    ctypeRNum  t256 = t122*t124*t215*(1.1E+1/1.2E+1);
    ctypeRNum  t257 = t123*t124*t215*(1.1E+1/1.2E+1);
    ctypeRNum  t258 = t101*t128*t201*(1.1E+1/2.4E+1);
    ctypeRNum  t259 = t102*t129*t201*(1.1E+1/2.4E+1);
    ctypeRNum  t260 = t103*t130*t201*(1.1E+1/2.4E+1);
    ctypeRNum  t261 = t104*t131*t203*(1.1E+1/2.4E+1);
    ctypeRNum  t262 = t105*t132*t203*(1.1E+1/2.4E+1);
    ctypeRNum  t263 = t106*t133*t203*(1.1E+1/2.4E+1);
    ctypeRNum  t264 = t107*t134*t205*(1.1E+1/2.4E+1);
    ctypeRNum  t265 = t108*t135*t205*(1.1E+1/2.4E+1);
    ctypeRNum  t266 = t109*t136*t205*(1.1E+1/2.4E+1);
    ctypeRNum  t267 = t110*t137*t207*(1.1E+1/2.4E+1);
    ctypeRNum  t268 = t111*t138*t207*(1.1E+1/2.4E+1);
    ctypeRNum  t269 = t112*t139*t207*(1.1E+1/2.4E+1);
    ctypeRNum  t270 = t113*t140*t209*(1.1E+1/2.4E+1);
    ctypeRNum  t271 = t114*t141*t209*(1.1E+1/2.4E+1);
    ctypeRNum  t272 = t115*t142*t209*(1.1E+1/2.4E+1);
    ctypeRNum  t273 = t116*t143*t211*(1.1E+1/2.4E+1);
    ctypeRNum  t274 = t117*t144*t211*(1.1E+1/2.4E+1);
    ctypeRNum  t275 = t118*t145*t211*(1.1E+1/2.4E+1);
    ctypeRNum  t276 = t119*t146*t213*(1.1E+1/2.4E+1);
    ctypeRNum  t277 = t120*t147*t213*(1.1E+1/2.4E+1);
    ctypeRNum  t278 = t121*t148*t213*(1.1E+1/2.4E+1);
    ctypeRNum  t279 = t122*t149*t215*(1.1E+1/2.4E+1);
    ctypeRNum  t280 = t123*t150*t215*(1.1E+1/2.4E+1);
    ctypeRNum  t281 = t124*t151*t215*(1.1E+1/2.4E+1);
    ctypeRNum  t282 = t187+t234;
    ctypeRNum  t283 = t188+t235;
    ctypeRNum  t284 = t189+t236;
    ctypeRNum  t285 = t226+t258+5.0E+1/3.0;
    ctypeRNum  t286 = t226+t259+5.0E+1/3.0;
    ctypeRNum  t287 = t226+t260+5.0E+1/3.0;
    ctypeRNum  t288 = t227+t261+5.0E+1/3.0;
    ctypeRNum  t289 = t227+t262+5.0E+1/3.0;
    ctypeRNum  t290 = t227+t263+5.0E+1/3.0;
    ctypeRNum  t291 = t228+t264+5.0E+1/3.0;
    ctypeRNum  t292 = t228+t265+5.0E+1/3.0;
    ctypeRNum  t293 = t228+t266+5.0E+1/3.0;
    ctypeRNum  t294 = t229+t267+5.0E+1/3.0;
    ctypeRNum  t295 = t229+t268+5.0E+1/3.0;
    ctypeRNum  t296 = t229+t269+5.0E+1/3.0;
    ctypeRNum  t297 = t230+t270+5.0E+1/3.0;
    ctypeRNum  t298 = t230+t271+5.0E+1/3.0;
    ctypeRNum  t299 = t230+t272+5.0E+1/3.0;
    ctypeRNum  t300 = t231+t273+5.0E+1/3.0;
    ctypeRNum  t301 = t231+t274+5.0E+1/3.0;
    ctypeRNum  t302 = t231+t275+5.0E+1/3.0;
    ctypeRNum  t303 = t232+t276+5.0E+1/3.0;
    ctypeRNum  t304 = t232+t277+5.0E+1/3.0;
    ctypeRNum  t305 = t232+t278+5.0E+1/3.0;
    ctypeRNum  t306 = t233+t279+5.0E+1/3.0;
    ctypeRNum  t307 = t233+t280+5.0E+1/3.0;
    ctypeRNum  t308 = t233+t281+5.0E+1/3.0;
    out[0] = t234*vec[10]+t235*vec[11]-t282*vec[4]-t283*vec[5]+t285*vec[9]-vec[3]*(t186+t226+t258+t32*t184*(1.1E+1/1.2E+1)+1.0E+2/3.0);
    out[1] = t234*vec[9]+t236*vec[11]-t282*vec[3]-t284*vec[5]+t286*vec[10]-vec[4]*(t186+t226+t259+t33*t184*(1.1E+1/1.2E+1)+1.0E+2/3.0);
    out[2] = t235*vec[9]+t236*vec[10]-t283*vec[3]-t284*vec[4]+t287*vec[11]-vec[5]*(t186+t226+t260+t34*t184*(1.1E+1/1.2E+1)+1.0E+2/3.0);
    out[3] = vec[0];
    out[4] = vec[1];
    out[5] = vec[2];
    out[6] = vec[9]*(-1.0E+2/3.0)+t218*vec[9]+t219*vec[9]+t234*vec[4]+t235*vec[5]+t237*vec[16]+t238*vec[17]+t285*vec[3]+t288*vec[15]-vec[10]*(t234+t105*t131*t203*(1.1E+1/2.4E+1))-vec[11]*(t235+t106*t131*t203*(1.1E+1/2.4E+1))-t101*t128*t201*vec[9]*(1.1E+1/2.4E+1)-t104*t131*t203*vec[9]*(1.1E+1/2.4E+1);
    out[7] = vec[10]*(-1.0E+2/3.0)+t218*vec[10]+t219*vec[10]+t234*vec[3]+t236*vec[5]+t237*vec[15]+t239*vec[17]+t286*vec[4]+t289*vec[16]-vec[9]*(t234+t104*t132*t203*(1.1E+1/2.4E+1))-vec[11]*(t236+t106*t132*t203*(1.1E+1/2.4E+1))-t102*t129*t201*vec[10]*(1.1E+1/2.4E+1)-t105*t132*t203*vec[10]*(1.1E+1/2.4E+1);
    out[8] = vec[11]*(-1.0E+2/3.0)+t218*vec[11]+t219*vec[11]+t235*vec[3]+t236*vec[4]+t238*vec[15]+t239*vec[16]+t287*vec[5]+t290*vec[17]-vec[9]*(t235+t104*t133*t203*(1.1E+1/2.4E+1))-vec[10]*(t236+t105*t133*t203*(1.1E+1/2.4E+1))-t103*t130*t201*vec[11]*(1.1E+1/2.4E+1)-t106*t133*t203*vec[11]*(1.1E+1/2.4E+1);
    out[9] = vec[6];
    out[10] = vec[7];
    out[11] = vec[8];
    out[12] = vec[15]*(-1.0E+2/3.0)+t219*vec[15]+t220*vec[15]+t237*vec[10]+t238*vec[11]+t240*vec[22]+t241*vec[23]+t288*vec[9]+t291*vec[21]-vec[16]*(t237+t108*t134*t205*(1.1E+1/2.4E+1))-vec[17]*(t238+t109*t134*t205*(1.1E+1/2.4E+1))-t104*t131*t203*vec[15]*(1.1E+1/2.4E+1)-t107*t134*t205*vec[15]*(1.1E+1/2.4E+1);
    out[13] = vec[16]*(-1.0E+2/3.0)+t219*vec[16]+t220*vec[16]+t237*vec[9]+t239*vec[11]+t240*vec[21]+t242*vec[23]+t289*vec[10]+t292*vec[22]-vec[15]*(t237+t107*t135*t205*(1.1E+1/2.4E+1))-vec[17]*(t239+t109*t135*t205*(1.1E+1/2.4E+1))-t105*t132*t203*vec[16]*(1.1E+1/2.4E+1)-t108*t135*t205*vec[16]*(1.1E+1/2.4E+1);
    out[14] = vec[17]*(-1.0E+2/3.0)+t219*vec[17]+t220*vec[17]+t238*vec[9]+t239*vec[10]+t241*vec[21]+t242*vec[22]+t290*vec[11]+t293*vec[23]-vec[15]*(t238+t107*t136*t205*(1.1E+1/2.4E+1))-vec[16]*(t239+t108*t136*t205*(1.1E+1/2.4E+1))-t106*t133*t203*vec[17]*(1.1E+1/2.4E+1)-t109*t136*t205*vec[17]*(1.1E+1/2.4E+1);
    out[15] = vec[12];
    out[16] = vec[13];
    out[17] = vec[14];
    out[18] = vec[21]*(-1.0E+2/3.0)+t220*vec[21]+t221*vec[21]+t240*vec[16]+t241*vec[17]+t243*vec[28]+t244*vec[29]+t291*vec[15]+t294*vec[27]-vec[22]*(t240+t111*t137*t207*(1.1E+1/2.4E+1))-vec[23]*(t241+t112*t137*t207*(1.1E+1/2.4E+1))-t107*t134*t205*vec[21]*(1.1E+1/2.4E+1)-t110*t137*t207*vec[21]*(1.1E+1/2.4E+1);
    out[19] = vec[22]*(-1.0E+2/3.0)+t220*vec[22]+t221*vec[22]+t240*vec[15]+t242*vec[17]+t243*vec[27]+t245*vec[29]+t292*vec[16]+t295*vec[28]-vec[21]*(t240+t110*t138*t207*(1.1E+1/2.4E+1))-vec[23]*(t242+t112*t138*t207*(1.1E+1/2.4E+1))-t108*t135*t205*vec[22]*(1.1E+1/2.4E+1)-t111*t138*t207*vec[22]*(1.1E+1/2.4E+1);
    out[20] = vec[23]*(-1.0E+2/3.0)+t220*vec[23]+t221*vec[23]+t241*vec[15]+t242*vec[16]+t244*vec[27]+t245*vec[28]+t293*vec[17]+t296*vec[29]-vec[21]*(t241+t110*t139*t207*(1.1E+1/2.4E+1))-vec[22]*(t242+t111*t139*t207*(1.1E+1/2.4E+1))-t109*t136*t205*vec[23]*(1.1E+1/2.4E+1)-t112*t139*t207*vec[23]*(1.1E+1/2.4E+1);
    out[21] = vec[18];
    out[22] = vec[19];
    out[23] = vec[20];
    out[24] = vec[27]*(-1.0E+2/3.0)+t221*vec[27]+t222*vec[27]+t243*vec[22]+t244*vec[23]+t246*vec[34]+t247*vec[35]+t294*vec[21]+t297*vec[33]-vec[28]*(t243+t114*t140*t209*(1.1E+1/2.4E+1))-vec[29]*(t244+t115*t140*t209*(1.1E+1/2.4E+1))-t110*t137*t207*vec[27]*(1.1E+1/2.4E+1)-t113*t140*t209*vec[27]*(1.1E+1/2.4E+1);
    out[25] = vec[28]*(-1.0E+2/3.0)+t221*vec[28]+t222*vec[28]+t243*vec[21]+t245*vec[23]+t246*vec[33]+t248*vec[35]+t295*vec[22]+t298*vec[34]-vec[27]*(t243+t113*t141*t209*(1.1E+1/2.4E+1))-vec[29]*(t245+t115*t141*t209*(1.1E+1/2.4E+1))-t111*t138*t207*vec[28]*(1.1E+1/2.4E+1)-t114*t141*t209*vec[28]*(1.1E+1/2.4E+1);
    out[26] = vec[29]*(-1.0E+2/3.0)+t221*vec[29]+t222*vec[29]+t244*vec[21]+t245*vec[22]+t247*vec[33]+t248*vec[34]+t296*vec[23]+t299*vec[35]-vec[27]*(t244+t113*t142*t209*(1.1E+1/2.4E+1))-vec[28]*(t245+t114*t142*t209*(1.1E+1/2.4E+1))-t112*t139*t207*vec[29]*(1.1E+1/2.4E+1)-t115*t142*t209*vec[29]*(1.1E+1/2.4E+1);
    out[27] = vec[24];
    out[28] = vec[25];
    out[29] = vec[26];
    out[30] = vec[33]*(-1.0E+2/3.0)+t222*vec[33]+t223*vec[33]+t246*vec[28]+t247*vec[29]+t249*vec[40]+t250*vec[41]+t297*vec[27]+t300*vec[39]-vec[34]*(t246+t117*t143*t211*(1.1E+1/2.4E+1))-vec[35]*(t247+t118*t143*t211*(1.1E+1/2.4E+1))-t113*t140*t209*vec[33]*(1.1E+1/2.4E+1)-t116*t143*t211*vec[33]*(1.1E+1/2.4E+1);
    out[31] = vec[34]*(-1.0E+2/3.0)+t222*vec[34]+t223*vec[34]+t246*vec[27]+t248*vec[29]+t249*vec[39]+t251*vec[41]+t298*vec[28]+t301*vec[40]-vec[33]*(t246+t116*t144*t211*(1.1E+1/2.4E+1))-vec[35]*(t248+t118*t144*t211*(1.1E+1/2.4E+1))-t114*t141*t209*vec[34]*(1.1E+1/2.4E+1)-t117*t144*t211*vec[34]*(1.1E+1/2.4E+1);
    out[32] = vec[35]*(-1.0E+2/3.0)+t222*vec[35]+t223*vec[35]+t247*vec[27]+t248*vec[28]+t250*vec[39]+t251*vec[40]+t299*vec[29]+t302*vec[41]-vec[33]*(t247+t116*t145*t211*(1.1E+1/2.4E+1))-vec[34]*(t248+t117*t145*t211*(1.1E+1/2.4E+1))-t115*t142*t209*vec[35]*(1.1E+1/2.4E+1)-t118*t145*t211*vec[35]*(1.1E+1/2.4E+1);
    out[33] = vec[30];
    out[34] = vec[31];
    out[35] = vec[32];
    out[36] = vec[39]*(-1.0E+2/3.0)+t223*vec[39]+t224*vec[39]+t249*vec[34]+t250*vec[35]+t252*vec[46]+t253*vec[47]+t300*vec[33]+t303*vec[45]-vec[40]*(t249+t120*t146*t213*(1.1E+1/2.4E+1))-vec[41]*(t250+t121*t146*t213*(1.1E+1/2.4E+1))-t116*t143*t211*vec[39]*(1.1E+1/2.4E+1)-t119*t146*t213*vec[39]*(1.1E+1/2.4E+1);
    out[37] = vec[40]*(-1.0E+2/3.0)+t223*vec[40]+t224*vec[40]+t249*vec[33]+t251*vec[35]+t252*vec[45]+t254*vec[47]+t301*vec[34]+t304*vec[46]-vec[39]*(t249+t119*t147*t213*(1.1E+1/2.4E+1))-vec[41]*(t251+t121*t147*t213*(1.1E+1/2.4E+1))-t117*t144*t211*vec[40]*(1.1E+1/2.4E+1)-t120*t147*t213*vec[40]*(1.1E+1/2.4E+1);
    out[38] = vec[41]*(-1.0E+2/3.0)+t223*vec[41]+t224*vec[41]+t250*vec[33]+t251*vec[34]+t253*vec[45]+t254*vec[46]+t302*vec[35]+t305*vec[47]-vec[39]*(t250+t119*t148*t213*(1.1E+1/2.4E+1))-vec[40]*(t251+t120*t148*t213*(1.1E+1/2.4E+1))-t118*t145*t211*vec[41]*(1.1E+1/2.4E+1)-t121*t148*t213*vec[41]*(1.1E+1/2.4E+1);
    out[39] = vec[36];
    out[40] = vec[37];
    out[41] = vec[38];
    out[42] = vec[45]*(-1.0E+2/3.0)+t224*vec[45]+t225*vec[45]+t252*vec[40]+t253*vec[41]+t255*vec[52]+t256*vec[53]+t303*vec[39]+t306*vec[51]-vec[46]*(t252+t123*t149*t215*(1.1E+1/2.4E+1))-vec[47]*(t253+t124*t149*t215*(1.1E+1/2.4E+1))-t119*t146*t213*vec[45]*(1.1E+1/2.4E+1)-t122*t149*t215*vec[45]*(1.1E+1/2.4E+1);
    out[43] = vec[46]*(-1.0E+2/3.0)+t224*vec[46]+t225*vec[46]+t252*vec[39]+t254*vec[41]+t255*vec[51]+t257*vec[53]+t304*vec[40]+t307*vec[52]-vec[45]*(t252+t122*t150*t215*(1.1E+1/2.4E+1))-vec[47]*(t254+t124*t150*t215*(1.1E+1/2.4E+1))-t120*t147*t213*vec[46]*(1.1E+1/2.4E+1)-t123*t150*t215*vec[46]*(1.1E+1/2.4E+1);
    out[44] = vec[47]*(-1.0E+2/3.0)+t224*vec[47]+t225*vec[47]+t253*vec[39]+t254*vec[40]+t256*vec[51]+t257*vec[52]+t305*vec[41]+t308*vec[53]-vec[45]*(t253+t122*t151*t215*(1.1E+1/2.4E+1))-vec[46]*(t254+t123*t151*t215*(1.1E+1/2.4E+1))-t121*t148*t213*vec[47]*(1.1E+1/2.4E+1)-t124*t151*t215*vec[47]*(1.1E+1/2.4E+1);
    out[45] = vec[42];
    out[46] = vec[43];
    out[47] = vec[44];
    out[48] = vec[51]*(-1.0E+2/3.0)+t216*vec[51]*(1.1E+1/1.2E+1)+t225*vec[51]+t255*vec[46]+t256*vec[47]+t306*vec[45]-vec[52]*(t255+t126*t152*t217*(1.1E+1/2.4E+1))-vec[53]*(t256+t127*t152*t217*(1.1E+1/2.4E+1))-t122*t149*t215*vec[51]*(1.1E+1/2.4E+1)-t125*t152*t217*vec[51]*(1.1E+1/2.4E+1);
    out[49] = vec[52]*(-1.0E+2/3.0)+t216*vec[52]*(1.1E+1/1.2E+1)+t225*vec[52]+t255*vec[45]+t257*vec[47]+t307*vec[46]-vec[51]*(t255+t125*t153*t217*(1.1E+1/2.4E+1))-vec[53]*(t257+t127*t153*t217*(1.1E+1/2.4E+1))-t123*t150*t215*vec[52]*(1.1E+1/2.4E+1)-t126*t153*t217*vec[52]*(1.1E+1/2.4E+1);
    out[50] = vec[53]*(-1.0E+2/3.0)+t216*vec[53]*(1.1E+1/1.2E+1)+t225*vec[53]+t256*vec[45]+t257*vec[46]+t308*vec[47]-vec[51]*(t256+t125*t154*t217*(1.1E+1/2.4E+1))-vec[52]*(t257+t126*t154*t217*(1.1E+1/2.4E+1))-t124*t151*t215*vec[53]*(1.1E+1/2.4E+1)-t127*t154*t217*vec[53]*(1.1E+1/2.4E+1);
    out[51] = vec[48];
    out[52] = vec[49];
    out[53] = vec[50];
    out[54] = (t217*(t35*vec[51]*1.1E+1+t199*vec[51]*2.0E+2+t93*x[48]+t94*x[48]-vec[51]*(t35+t180+t181)*1.1E+1-vec[52]*x[48]*x[55]*1.1E+1-vec[53]*x[48]*x[56]*1.1E+1))/1.2E+1-(t217*x[54]*(t93+t94+t99+t100))/1.2E+1;
    out[55] = (t217*(t36*vec[52]*1.1E+1+t199*vec[52]*2.0E+2+t92*x[49]+t94*x[49]-vec[52]*(t36+t179+t181)*1.1E+1-vec[51]*x[49]*x[54]*1.1E+1-vec[53]*x[49]*x[56]*1.1E+1))/1.2E+1-(t217*x[55]*(t92+t94+t98+t100))/1.2E+1;
    out[56] = (t217*(t37*vec[53]*1.1E+1+t199*vec[53]*2.0E+2+t92*x[50]+t93*x[50]-vec[53]*(t37+t179+t180)*1.1E+1-vec[51]*x[50]*x[54]*1.1E+1-vec[52]*x[50]*x[55]*1.1E+1))/1.2E+1-(t217*x[56]*(t92+t93+t98+t99))/1.2E+1;
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Chain_10_ProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = vec[54];
    out[1] = vec[55];
    out[2] = vec[56];
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void Chain_10_ProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void Chain_10_ProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*(u[0]*u[0])+pCost_[2]*(u[1]*u[1])+pCost_[2]*(u[2]*u[2])+pCost_[1]*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+pCost_[1]*(x[9]*x[9]+x[10]*x[10]+x[11]*x[11])+pCost_[1]*(x[15]*x[15]+x[16]*x[16]+x[17]*x[17])+pCost_[1]*(x[21]*x[21]+x[22]*x[22]+x[23]*x[23])+pCost_[1]*(x[27]*x[27]+x[28]*x[28]+x[29]*x[29])+pCost_[1]*(x[33]*x[33]+x[34]*x[34]+x[35]*x[35])+pCost_[1]*(x[39]*x[39]+x[40]*x[40]+x[41]*x[41])+pCost_[1]*(x[45]*x[45]+x[46]*x[46]+x[47]*x[47])+pCost_[1]*(x[51]*x[51]+x[52]*x[52]+x[53]*x[53])+pCost_[0]*(POW2(x[54]-5.0)+x[55]*x[55]+x[56]*x[56]);
}


/** Gradient dl/dx **/
void Chain_10_ProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
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
    out[54] = pCost_[0]*(x[54]-5.0)*2.0;
    out[55] = pCost_[0]*x[55]*2.0;
    out[56] = pCost_[0]*x[56]*2.0;
}


/** Gradient dl/du **/
void Chain_10_ProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*u[0]*2.0;
    out[1] = pCost_[2]*u[1]*2.0;
    out[2] = pCost_[2]*u[2]*2.0;
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Chain_10_ProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Gradient dV/dx **/
void Chain_10_ProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
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
}


/** Gradient dV/dT **/
void Chain_10_ProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void Chain_10_ProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Chain_10_ProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Chain_10_ProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}



/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void Chain_10_ProblemDescription::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void Chain_10_ProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void Chain_10_ProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}