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
* This probfct-file describes the nonlinear chain problem with 14 chain elements from
* Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
* Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
* Process Control 22(3), 540-550 (2012)
*
*/

#include "NLChain_14_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

using namespace grampc;


Chain_14_ProblemDescription::Chain_14_ProblemDescription(const std::vector<typeRNum>& pCost)
: pCost_(pCost)
{
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void Chain_14_ProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 81;
    *Nu  = 3;
    *Np  = 0;
    *Nh  = 0;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void Chain_14_ProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
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
    ctypeRNum  t38 = -x[72];
    ctypeRNum  t39 = -x[73];
    ctypeRNum  t40 = -x[74];
    ctypeRNum  t41 = -x[78];
    ctypeRNum  t42 = -x[79];
    ctypeRNum  t43 = -x[80];
    ctypeRNum  t44 = t5+x[0];
    ctypeRNum  t45 = t6+x[1];
    ctypeRNum  t46 = t7+x[2];
    ctypeRNum  t47 = t8+x[6];
    ctypeRNum  t48 = t9+x[7];
    ctypeRNum  t49 = t10+x[8];
    ctypeRNum  t50 = t11+x[12];
    ctypeRNum  t51 = t12+x[13];
    ctypeRNum  t52 = t13+x[14];
    ctypeRNum  t53 = t14+x[18];
    ctypeRNum  t54 = t15+x[19];
    ctypeRNum  t55 = t16+x[20];
    ctypeRNum  t56 = t17+x[24];
    ctypeRNum  t57 = t18+x[25];
    ctypeRNum  t58 = t19+x[26];
    ctypeRNum  t59 = t20+x[30];
    ctypeRNum  t60 = t21+x[31];
    ctypeRNum  t61 = t22+x[32];
    ctypeRNum  t62 = t23+x[36];
    ctypeRNum  t63 = t24+x[37];
    ctypeRNum  t64 = t25+x[38];
    ctypeRNum  t65 = t26+x[42];
    ctypeRNum  t66 = t27+x[43];
    ctypeRNum  t67 = t28+x[44];
    ctypeRNum  t68 = t29+x[48];
    ctypeRNum  t69 = t30+x[49];
    ctypeRNum  t70 = t31+x[50];
    ctypeRNum  t71 = t32+x[54];
    ctypeRNum  t72 = t33+x[55];
    ctypeRNum  t73 = t34+x[56];
    ctypeRNum  t74 = t35+x[60];
    ctypeRNum  t75 = t36+x[61];
    ctypeRNum  t76 = t37+x[62];
    ctypeRNum  t77 = t38+x[66];
    ctypeRNum  t78 = t39+x[67];
    ctypeRNum  t79 = t40+x[68];
    ctypeRNum  t80 = t41+x[72];
    ctypeRNum  t81 = t42+x[73];
    ctypeRNum  t82 = t43+x[74];
    ctypeRNum  t122 = t2+t3+t4;
    ctypeRNum  t83 = t44*t44;
    ctypeRNum  t84 = t45*t45;
    ctypeRNum  t85 = t46*t46;
    ctypeRNum  t86 = t47*t47;
    ctypeRNum  t87 = t48*t48;
    ctypeRNum  t88 = t49*t49;
    ctypeRNum  t89 = t50*t50;
    ctypeRNum  t90 = t51*t51;
    ctypeRNum  t91 = t52*t52;
    ctypeRNum  t92 = t53*t53;
    ctypeRNum  t93 = t54*t54;
    ctypeRNum  t94 = t55*t55;
    ctypeRNum  t95 = t56*t56;
    ctypeRNum  t96 = t57*t57;
    ctypeRNum  t97 = t58*t58;
    ctypeRNum  t98 = t59*t59;
    ctypeRNum  t99 = t60*t60;
    ctypeRNum  t100 = t61*t61;
    ctypeRNum  t101 = t62*t62;
    ctypeRNum  t102 = t63*t63;
    ctypeRNum  t103 = t64*t64;
    ctypeRNum  t104 = t65*t65;
    ctypeRNum  t105 = t66*t66;
    ctypeRNum  t106 = t67*t67;
    ctypeRNum  t107 = t68*t68;
    ctypeRNum  t108 = t69*t69;
    ctypeRNum  t109 = t70*t70;
    ctypeRNum  t110 = t71*t71;
    ctypeRNum  t111 = t72*t72;
    ctypeRNum  t112 = t73*t73;
    ctypeRNum  t113 = t74*t74;
    ctypeRNum  t114 = t75*t75;
    ctypeRNum  t115 = t76*t76;
    ctypeRNum  t116 = t77*t77;
    ctypeRNum  t117 = t78*t78;
    ctypeRNum  t118 = t79*t79;
    ctypeRNum  t119 = t80*t80;
    ctypeRNum  t120 = t81*t81;
    ctypeRNum  t121 = t82*t82;
    ctypeRNum  t123 = 1.0/sqrt(t122);
    ctypeRNum  t124 = t123*(1.1E+1/4.0E+2);
    ctypeRNum  t125 = t83+t84+t85;
    ctypeRNum  t126 = t86+t87+t88;
    ctypeRNum  t127 = t89+t90+t91;
    ctypeRNum  t128 = t92+t93+t94;
    ctypeRNum  t129 = t95+t96+t97;
    ctypeRNum  t130 = t98+t99+t100;
    ctypeRNum  t131 = t101+t102+t103;
    ctypeRNum  t132 = t104+t105+t106;
    ctypeRNum  t133 = t107+t108+t109;
    ctypeRNum  t134 = t110+t111+t112;
    ctypeRNum  t135 = t113+t114+t115;
    ctypeRNum  t136 = t116+t117+t118;
    ctypeRNum  t137 = t119+t120+t121;
    ctypeRNum  t138 = t124-7.0/1.0E+1;
    ctypeRNum  t139 = 1.0/sqrt(t125);
    ctypeRNum  t140 = 1.0/sqrt(t126);
    ctypeRNum  t141 = 1.0/sqrt(t127);
    ctypeRNum  t142 = 1.0/sqrt(t128);
    ctypeRNum  t143 = 1.0/sqrt(t129);
    ctypeRNum  t144 = 1.0/sqrt(t130);
    ctypeRNum  t145 = 1.0/sqrt(t131);
    ctypeRNum  t146 = 1.0/sqrt(t132);
    ctypeRNum  t147 = 1.0/sqrt(t133);
    ctypeRNum  t148 = 1.0/sqrt(t134);
    ctypeRNum  t149 = 1.0/sqrt(t135);
    ctypeRNum  t150 = 1.0/sqrt(t136);
    ctypeRNum  t151 = 1.0/sqrt(t137);
    ctypeRNum  t152 = t139*(1.1E+1/4.0E+2);
    ctypeRNum  t153 = t140*(1.1E+1/4.0E+2);
    ctypeRNum  t154 = t141*(1.1E+1/4.0E+2);
    ctypeRNum  t155 = t142*(1.1E+1/4.0E+2);
    ctypeRNum  t156 = t143*(1.1E+1/4.0E+2);
    ctypeRNum  t157 = t144*(1.1E+1/4.0E+2);
    ctypeRNum  t158 = t145*(1.1E+1/4.0E+2);
    ctypeRNum  t159 = t146*(1.1E+1/4.0E+2);
    ctypeRNum  t160 = t147*(1.1E+1/4.0E+2);
    ctypeRNum  t161 = t148*(1.1E+1/4.0E+2);
    ctypeRNum  t162 = t149*(1.1E+1/4.0E+2);
    ctypeRNum  t163 = t150*(1.1E+1/4.0E+2);
    ctypeRNum  t164 = t151*(1.1E+1/4.0E+2);
    ctypeRNum  t165 = t152-7.0/1.0E+1;
    ctypeRNum  t166 = t153-7.0/1.0E+1;
    ctypeRNum  t167 = t154-7.0/1.0E+1;
    ctypeRNum  t168 = t155-7.0/1.0E+1;
    ctypeRNum  t169 = t156-7.0/1.0E+1;
    ctypeRNum  t170 = t157-7.0/1.0E+1;
    ctypeRNum  t171 = t158-7.0/1.0E+1;
    ctypeRNum  t172 = t159-7.0/1.0E+1;
    ctypeRNum  t173 = t160-7.0/1.0E+1;
    ctypeRNum  t174 = t161-7.0/1.0E+1;
    ctypeRNum  t175 = t162-7.0/1.0E+1;
    ctypeRNum  t176 = t163-7.0/1.0E+1;
    ctypeRNum  t177 = t164-7.0/1.0E+1;
    ctypeRNum  t178 = t44*t165*(1.0E+2/3.0);
    ctypeRNum  t179 = t45*t165*(1.0E+2/3.0);
    ctypeRNum  t180 = t46*t165*(1.0E+2/3.0);
    ctypeRNum  t181 = t47*t166*(1.0E+2/3.0);
    ctypeRNum  t182 = t48*t166*(1.0E+2/3.0);
    ctypeRNum  t183 = t49*t166*(1.0E+2/3.0);
    ctypeRNum  t184 = t50*t167*(1.0E+2/3.0);
    ctypeRNum  t185 = t51*t167*(1.0E+2/3.0);
    ctypeRNum  t186 = t52*t167*(1.0E+2/3.0);
    ctypeRNum  t187 = t53*t168*(1.0E+2/3.0);
    ctypeRNum  t188 = t54*t168*(1.0E+2/3.0);
    ctypeRNum  t189 = t55*t168*(1.0E+2/3.0);
    ctypeRNum  t190 = t56*t169*(1.0E+2/3.0);
    ctypeRNum  t191 = t57*t169*(1.0E+2/3.0);
    ctypeRNum  t192 = t58*t169*(1.0E+2/3.0);
    ctypeRNum  t193 = t59*t170*(1.0E+2/3.0);
    ctypeRNum  t194 = t60*t170*(1.0E+2/3.0);
    ctypeRNum  t195 = t61*t170*(1.0E+2/3.0);
    ctypeRNum  t196 = t62*t171*(1.0E+2/3.0);
    ctypeRNum  t197 = t63*t171*(1.0E+2/3.0);
    ctypeRNum  t198 = t64*t171*(1.0E+2/3.0);
    ctypeRNum  t199 = t65*t172*(1.0E+2/3.0);
    ctypeRNum  t200 = t66*t172*(1.0E+2/3.0);
    ctypeRNum  t201 = t67*t172*(1.0E+2/3.0);
    ctypeRNum  t202 = t68*t173*(1.0E+2/3.0);
    ctypeRNum  t203 = t69*t173*(1.0E+2/3.0);
    ctypeRNum  t204 = t70*t173*(1.0E+2/3.0);
    ctypeRNum  t205 = t71*t174*(1.0E+2/3.0);
    ctypeRNum  t206 = t72*t174*(1.0E+2/3.0);
    ctypeRNum  t207 = t73*t174*(1.0E+2/3.0);
    ctypeRNum  t208 = t74*t175*(1.0E+2/3.0);
    ctypeRNum  t209 = t75*t175*(1.0E+2/3.0);
    ctypeRNum  t210 = t76*t175*(1.0E+2/3.0);
    ctypeRNum  t211 = t77*t176*(1.0E+2/3.0);
    ctypeRNum  t212 = t78*t176*(1.0E+2/3.0);
    ctypeRNum  t213 = t79*t176*(1.0E+2/3.0);
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = t178+t138*x[0]*(1.0E+2/3.0);
    out[4] = t179+t138*x[1]*(1.0E+2/3.0);
    out[5] = t180+t138*x[2]*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[6] = x[9];
    out[7] = x[10];
    out[8] = x[11];
    out[9] = -t178+t181;
    out[10] = -t179+t182;
    out[11] = -t180+t183-9.81E+2/1.0E+2;
    out[12] = x[15];
    out[13] = x[16];
    out[14] = x[17];
    out[15] = -t181+t184;
    out[16] = -t182+t185;
    out[17] = -t183+t186-9.81E+2/1.0E+2;
    out[18] = x[21];
    out[19] = x[22];
    out[20] = x[23];
    out[21] = -t184+t187;
    out[22] = -t185+t188;
    out[23] = -t186+t189-9.81E+2/1.0E+2;
    out[24] = x[27];
    out[25] = x[28];
    out[26] = x[29];
    out[27] = -t187+t190;
    out[28] = -t188+t191;
    out[29] = -t189+t192-9.81E+2/1.0E+2;
    out[30] = x[33];
    out[31] = x[34];
    out[32] = x[35];
    out[33] = -t190+t193;
    out[34] = -t191+t194;
    out[35] = -t192+t195-9.81E+2/1.0E+2;
    out[36] = x[39];
    out[37] = x[40];
    out[38] = x[41];
    out[39] = -t193+t196;
    out[40] = -t194+t197;
    out[41] = -t195+t198-9.81E+2/1.0E+2;
    out[42] = x[45];
    out[43] = x[46];
    out[44] = x[47];
    out[45] = -t196+t199;
    out[46] = -t197+t200;
    out[47] = -t198+t201-9.81E+2/1.0E+2;
    out[48] = x[51];
    out[49] = x[52];
    out[50] = x[53];
    out[51] = -t199+t202;
    out[52] = -t200+t203;
    out[53] = -t201+t204-9.81E+2/1.0E+2;
    out[54] = x[57];
    out[55] = x[58];
    out[56] = x[59];
    out[57] = -t202+t205;
    out[58] = -t203+t206;
    out[59] = -t204+t207-9.81E+2/1.0E+2;
    out[60] = x[63];
    out[61] = x[64];
    out[62] = x[65];
    out[63] = -t205+t208;
    out[64] = -t206+t209;
    out[65] = -t207+t210-9.81E+2/1.0E+2;
    out[66] = x[69];
    out[67] = x[70];
    out[68] = x[71];
    out[69] = -t208+t211;
    out[70] = -t209+t212;
    out[71] = -t210+t213-9.81E+2/1.0E+2;
    out[72] = x[75];
    out[73] = x[76];
    out[74] = x[77];
    out[75] = -t211+t80*t177*(1.0E+2/3.0);
    out[76] = -t212+t81*t177*(1.0E+2/3.0);
    out[77] = -t213+t82*t177*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[78] = u[0];
    out[79] = u[1];
    out[80] = u[2];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Chain_14_ProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
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
    ctypeRNum  t38 = x[72]*2.0;
    ctypeRNum  t39 = x[73]*2.0;
    ctypeRNum  t40 = x[74]*2.0;
    ctypeRNum  t41 = x[78]*2.0;
    ctypeRNum  t42 = x[79]*2.0;
    ctypeRNum  t43 = x[80]*2.0;
    ctypeRNum  t44 = x[0]*x[0];
    ctypeRNum  t45 = x[1]*x[1];
    ctypeRNum  t46 = x[2]*x[2];
    ctypeRNum  t47 = x[72]*x[72];
    ctypeRNum  t48 = x[73]*x[73];
    ctypeRNum  t49 = x[74]*x[74];
    ctypeRNum  t50 = -x[6];
    ctypeRNum  t52 = -x[7];
    ctypeRNum  t54 = -x[8];
    ctypeRNum  t56 = -x[12];
    ctypeRNum  t58 = -x[13];
    ctypeRNum  t60 = -x[14];
    ctypeRNum  t62 = -x[18];
    ctypeRNum  t64 = -x[19];
    ctypeRNum  t66 = -x[20];
    ctypeRNum  t68 = -x[24];
    ctypeRNum  t70 = -x[25];
    ctypeRNum  t72 = -x[26];
    ctypeRNum  t74 = -x[30];
    ctypeRNum  t76 = -x[31];
    ctypeRNum  t78 = -x[32];
    ctypeRNum  t80 = -x[36];
    ctypeRNum  t82 = -x[37];
    ctypeRNum  t84 = -x[38];
    ctypeRNum  t86 = -x[42];
    ctypeRNum  t88 = -x[43];
    ctypeRNum  t90 = -x[44];
    ctypeRNum  t92 = -x[48];
    ctypeRNum  t94 = -x[49];
    ctypeRNum  t96 = -x[50];
    ctypeRNum  t98 = -x[54];
    ctypeRNum  t100 = -x[55];
    ctypeRNum  t102 = -x[56];
    ctypeRNum  t104 = -x[60];
    ctypeRNum  t106 = -x[61];
    ctypeRNum  t108 = -x[62];
    ctypeRNum  t110 = -x[66];
    ctypeRNum  t112 = -x[67];
    ctypeRNum  t114 = -x[68];
    ctypeRNum  t116 = -x[72];
    ctypeRNum  t118 = -x[73];
    ctypeRNum  t120 = -x[74];
    ctypeRNum  t122 = -x[78];
    ctypeRNum  t124 = -x[79];
    ctypeRNum  t126 = -x[80];
    ctypeRNum  t128 = vec[75]*x[72]*1.1E+1;
    ctypeRNum  t129 = vec[76]*x[73]*1.1E+1;
    ctypeRNum  t130 = vec[77]*x[74]*1.1E+1;
    ctypeRNum  t131 = vec[75]*x[78]*1.1E+1;
    ctypeRNum  t132 = vec[76]*x[79]*1.1E+1;
    ctypeRNum  t133 = vec[77]*x[80]*1.1E+1;
    ctypeRNum  t51 = -t5;
    ctypeRNum  t53 = -t6;
    ctypeRNum  t55 = -t7;
    ctypeRNum  t57 = -t8;
    ctypeRNum  t59 = -t9;
    ctypeRNum  t61 = -t10;
    ctypeRNum  t63 = -t11;
    ctypeRNum  t65 = -t12;
    ctypeRNum  t67 = -t13;
    ctypeRNum  t69 = -t14;
    ctypeRNum  t71 = -t15;
    ctypeRNum  t73 = -t16;
    ctypeRNum  t75 = -t17;
    ctypeRNum  t77 = -t18;
    ctypeRNum  t79 = -t19;
    ctypeRNum  t81 = -t20;
    ctypeRNum  t83 = -t21;
    ctypeRNum  t85 = -t22;
    ctypeRNum  t87 = -t23;
    ctypeRNum  t89 = -t24;
    ctypeRNum  t91 = -t25;
    ctypeRNum  t93 = -t26;
    ctypeRNum  t95 = -t27;
    ctypeRNum  t97 = -t28;
    ctypeRNum  t99 = -t29;
    ctypeRNum  t101 = -t30;
    ctypeRNum  t103 = -t31;
    ctypeRNum  t105 = -t32;
    ctypeRNum  t107 = -t33;
    ctypeRNum  t109 = -t34;
    ctypeRNum  t111 = -t35;
    ctypeRNum  t113 = -t36;
    ctypeRNum  t115 = -t37;
    ctypeRNum  t117 = -t38;
    ctypeRNum  t119 = -t39;
    ctypeRNum  t121 = -t40;
    ctypeRNum  t123 = -t41;
    ctypeRNum  t125 = -t42;
    ctypeRNum  t127 = -t43;
    ctypeRNum  t134 = -t131;
    ctypeRNum  t135 = -t132;
    ctypeRNum  t136 = -t133;
    ctypeRNum  t137 = t50+x[0];
    ctypeRNum  t138 = t52+x[1];
    ctypeRNum  t139 = t54+x[2];
    ctypeRNum  t140 = t56+x[6];
    ctypeRNum  t141 = t58+x[7];
    ctypeRNum  t142 = t60+x[8];
    ctypeRNum  t143 = t62+x[12];
    ctypeRNum  t144 = t64+x[13];
    ctypeRNum  t145 = t66+x[14];
    ctypeRNum  t146 = t68+x[18];
    ctypeRNum  t147 = t70+x[19];
    ctypeRNum  t148 = t72+x[20];
    ctypeRNum  t149 = t74+x[24];
    ctypeRNum  t150 = t76+x[25];
    ctypeRNum  t151 = t78+x[26];
    ctypeRNum  t152 = t80+x[30];
    ctypeRNum  t153 = t82+x[31];
    ctypeRNum  t154 = t84+x[32];
    ctypeRNum  t155 = t86+x[36];
    ctypeRNum  t156 = t88+x[37];
    ctypeRNum  t157 = t90+x[38];
    ctypeRNum  t158 = t92+x[42];
    ctypeRNum  t159 = t94+x[43];
    ctypeRNum  t160 = t96+x[44];
    ctypeRNum  t161 = t98+x[48];
    ctypeRNum  t162 = t100+x[49];
    ctypeRNum  t163 = t102+x[50];
    ctypeRNum  t164 = t104+x[54];
    ctypeRNum  t165 = t106+x[55];
    ctypeRNum  t166 = t108+x[56];
    ctypeRNum  t167 = t110+x[60];
    ctypeRNum  t168 = t112+x[61];
    ctypeRNum  t169 = t114+x[62];
    ctypeRNum  t170 = t116+x[66];
    ctypeRNum  t171 = t118+x[67];
    ctypeRNum  t172 = t120+x[68];
    ctypeRNum  t173 = t122+x[72];
    ctypeRNum  t174 = t124+x[73];
    ctypeRNum  t175 = t126+x[74];
    ctypeRNum  t254 = t44+t45+t46;
    ctypeRNum  t176 = t2+t51;
    ctypeRNum  t177 = t3+t53;
    ctypeRNum  t178 = t4+t55;
    ctypeRNum  t179 = t5+t57;
    ctypeRNum  t180 = t6+t59;
    ctypeRNum  t181 = t7+t61;
    ctypeRNum  t182 = t8+t63;
    ctypeRNum  t183 = t9+t65;
    ctypeRNum  t184 = t10+t67;
    ctypeRNum  t185 = t11+t69;
    ctypeRNum  t186 = t12+t71;
    ctypeRNum  t187 = t13+t73;
    ctypeRNum  t188 = t14+t75;
    ctypeRNum  t189 = t15+t77;
    ctypeRNum  t190 = t16+t79;
    ctypeRNum  t191 = t17+t81;
    ctypeRNum  t192 = t18+t83;
    ctypeRNum  t193 = t19+t85;
    ctypeRNum  t194 = t20+t87;
    ctypeRNum  t195 = t21+t89;
    ctypeRNum  t196 = t22+t91;
    ctypeRNum  t197 = t23+t93;
    ctypeRNum  t198 = t24+t95;
    ctypeRNum  t199 = t25+t97;
    ctypeRNum  t200 = t26+t99;
    ctypeRNum  t201 = t27+t101;
    ctypeRNum  t202 = t28+t103;
    ctypeRNum  t203 = t29+t105;
    ctypeRNum  t204 = t30+t107;
    ctypeRNum  t205 = t31+t109;
    ctypeRNum  t206 = t32+t111;
    ctypeRNum  t207 = t33+t113;
    ctypeRNum  t208 = t34+t115;
    ctypeRNum  t209 = t35+t117;
    ctypeRNum  t210 = t36+t119;
    ctypeRNum  t211 = t37+t121;
    ctypeRNum  t212 = t38+t123;
    ctypeRNum  t213 = t39+t125;
    ctypeRNum  t214 = t40+t127;
    ctypeRNum  t215 = t137*t137;
    ctypeRNum  t216 = t138*t138;
    ctypeRNum  t217 = t139*t139;
    ctypeRNum  t218 = t140*t140;
    ctypeRNum  t219 = t141*t141;
    ctypeRNum  t220 = t142*t142;
    ctypeRNum  t221 = t143*t143;
    ctypeRNum  t222 = t144*t144;
    ctypeRNum  t223 = t145*t145;
    ctypeRNum  t224 = t146*t146;
    ctypeRNum  t225 = t147*t147;
    ctypeRNum  t226 = t148*t148;
    ctypeRNum  t227 = t149*t149;
    ctypeRNum  t228 = t150*t150;
    ctypeRNum  t229 = t151*t151;
    ctypeRNum  t230 = t152*t152;
    ctypeRNum  t231 = t153*t153;
    ctypeRNum  t232 = t154*t154;
    ctypeRNum  t233 = t155*t155;
    ctypeRNum  t234 = t156*t156;
    ctypeRNum  t235 = t157*t157;
    ctypeRNum  t236 = t158*t158;
    ctypeRNum  t237 = t159*t159;
    ctypeRNum  t238 = t160*t160;
    ctypeRNum  t239 = t161*t161;
    ctypeRNum  t240 = t162*t162;
    ctypeRNum  t241 = t163*t163;
    ctypeRNum  t242 = t164*t164;
    ctypeRNum  t243 = t165*t165;
    ctypeRNum  t244 = t166*t166;
    ctypeRNum  t245 = t167*t167;
    ctypeRNum  t246 = t168*t168;
    ctypeRNum  t247 = t169*t169;
    ctypeRNum  t248 = t170*t170;
    ctypeRNum  t249 = t171*t171;
    ctypeRNum  t250 = t172*t172;
    ctypeRNum  t251 = t173*t173;
    ctypeRNum  t252 = t174*t174;
    ctypeRNum  t253 = t175*t175;
    ctypeRNum  t255 = 1.0/sqrt(t254);
    ctypeRNum  t256 = t255*t255*t255;
    ctypeRNum  t257 = t255*(1.1E+1/1.2E+1);
    ctypeRNum  t262 = t215+t216+t217;
    ctypeRNum  t263 = t218+t219+t220;
    ctypeRNum  t264 = t221+t222+t223;
    ctypeRNum  t265 = t224+t225+t226;
    ctypeRNum  t266 = t227+t228+t229;
    ctypeRNum  t267 = t230+t231+t232;
    ctypeRNum  t268 = t233+t234+t235;
    ctypeRNum  t269 = t236+t237+t238;
    ctypeRNum  t270 = t239+t240+t241;
    ctypeRNum  t271 = t242+t243+t244;
    ctypeRNum  t272 = t245+t246+t247;
    ctypeRNum  t273 = t248+t249+t250;
    ctypeRNum  t274 = t251+t252+t253;
    ctypeRNum  t258 = -t257;
    ctypeRNum  t259 = t256*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t260 = t256*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t261 = t256*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t275 = pow(t274,3.0/2.0);
    ctypeRNum  t276 = 1.0/sqrt(t262);
    ctypeRNum  t278 = 1.0/sqrt(t263);
    ctypeRNum  t280 = 1.0/sqrt(t264);
    ctypeRNum  t282 = 1.0/sqrt(t265);
    ctypeRNum  t284 = 1.0/sqrt(t266);
    ctypeRNum  t286 = 1.0/sqrt(t267);
    ctypeRNum  t288 = 1.0/sqrt(t268);
    ctypeRNum  t290 = 1.0/sqrt(t269);
    ctypeRNum  t292 = 1.0/sqrt(t270);
    ctypeRNum  t294 = 1.0/sqrt(t271);
    ctypeRNum  t296 = 1.0/sqrt(t272);
    ctypeRNum  t298 = 1.0/sqrt(t273);
    ctypeRNum  t300 = 1.0/sqrt(t274);
    ctypeRNum  t277 = t276*t276*t276;
    ctypeRNum  t279 = t278*t278*t278;
    ctypeRNum  t281 = t280*t280*t280;
    ctypeRNum  t283 = t282*t282*t282;
    ctypeRNum  t285 = t284*t284*t284;
    ctypeRNum  t287 = t286*t286*t286;
    ctypeRNum  t289 = t288*t288*t288;
    ctypeRNum  t291 = t290*t290*t290;
    ctypeRNum  t293 = t292*t292*t292;
    ctypeRNum  t295 = t294*t294*t294;
    ctypeRNum  t297 = t296*t296*t296;
    ctypeRNum  t299 = t298*t298*t298;
    ctypeRNum  t301 = 1.0/t275;
    ctypeRNum  t302 = t276*(1.1E+1/1.2E+1);
    ctypeRNum  t303 = t278*(1.1E+1/1.2E+1);
    ctypeRNum  t304 = t280*(1.1E+1/1.2E+1);
    ctypeRNum  t305 = t282*(1.1E+1/1.2E+1);
    ctypeRNum  t306 = t284*(1.1E+1/1.2E+1);
    ctypeRNum  t307 = t286*(1.1E+1/1.2E+1);
    ctypeRNum  t308 = t288*(1.1E+1/1.2E+1);
    ctypeRNum  t309 = t290*(1.1E+1/1.2E+1);
    ctypeRNum  t310 = t292*(1.1E+1/1.2E+1);
    ctypeRNum  t311 = t294*(1.1E+1/1.2E+1);
    ctypeRNum  t312 = t296*(1.1E+1/1.2E+1);
    ctypeRNum  t313 = t298*(1.1E+1/1.2E+1);
    ctypeRNum  t314 = -t302;
    ctypeRNum  t315 = -t303;
    ctypeRNum  t316 = -t304;
    ctypeRNum  t317 = -t305;
    ctypeRNum  t318 = -t306;
    ctypeRNum  t319 = -t307;
    ctypeRNum  t320 = -t308;
    ctypeRNum  t321 = -t309;
    ctypeRNum  t322 = -t310;
    ctypeRNum  t323 = -t311;
    ctypeRNum  t324 = -t312;
    ctypeRNum  t325 = -t313;
    ctypeRNum  t326 = t137*t138*t277*(1.1E+1/1.2E+1);
    ctypeRNum  t327 = t137*t139*t277*(1.1E+1/1.2E+1);
    ctypeRNum  t328 = t138*t139*t277*(1.1E+1/1.2E+1);
    ctypeRNum  t329 = t140*t141*t279*(1.1E+1/1.2E+1);
    ctypeRNum  t330 = t140*t142*t279*(1.1E+1/1.2E+1);
    ctypeRNum  t331 = t141*t142*t279*(1.1E+1/1.2E+1);
    ctypeRNum  t332 = t143*t144*t281*(1.1E+1/1.2E+1);
    ctypeRNum  t333 = t143*t145*t281*(1.1E+1/1.2E+1);
    ctypeRNum  t334 = t144*t145*t281*(1.1E+1/1.2E+1);
    ctypeRNum  t335 = t146*t147*t283*(1.1E+1/1.2E+1);
    ctypeRNum  t336 = t146*t148*t283*(1.1E+1/1.2E+1);
    ctypeRNum  t337 = t147*t148*t283*(1.1E+1/1.2E+1);
    ctypeRNum  t338 = t149*t150*t285*(1.1E+1/1.2E+1);
    ctypeRNum  t339 = t149*t151*t285*(1.1E+1/1.2E+1);
    ctypeRNum  t340 = t150*t151*t285*(1.1E+1/1.2E+1);
    ctypeRNum  t341 = t152*t153*t287*(1.1E+1/1.2E+1);
    ctypeRNum  t342 = t152*t154*t287*(1.1E+1/1.2E+1);
    ctypeRNum  t343 = t153*t154*t287*(1.1E+1/1.2E+1);
    ctypeRNum  t344 = t155*t156*t289*(1.1E+1/1.2E+1);
    ctypeRNum  t345 = t155*t157*t289*(1.1E+1/1.2E+1);
    ctypeRNum  t346 = t156*t157*t289*(1.1E+1/1.2E+1);
    ctypeRNum  t347 = t158*t159*t291*(1.1E+1/1.2E+1);
    ctypeRNum  t348 = t158*t160*t291*(1.1E+1/1.2E+1);
    ctypeRNum  t349 = t159*t160*t291*(1.1E+1/1.2E+1);
    ctypeRNum  t350 = t161*t162*t293*(1.1E+1/1.2E+1);
    ctypeRNum  t351 = t161*t163*t293*(1.1E+1/1.2E+1);
    ctypeRNum  t352 = t162*t163*t293*(1.1E+1/1.2E+1);
    ctypeRNum  t353 = t164*t165*t295*(1.1E+1/1.2E+1);
    ctypeRNum  t354 = t164*t166*t295*(1.1E+1/1.2E+1);
    ctypeRNum  t355 = t165*t166*t295*(1.1E+1/1.2E+1);
    ctypeRNum  t356 = t167*t168*t297*(1.1E+1/1.2E+1);
    ctypeRNum  t357 = t167*t169*t297*(1.1E+1/1.2E+1);
    ctypeRNum  t358 = t168*t169*t297*(1.1E+1/1.2E+1);
    ctypeRNum  t359 = t170*t171*t299*(1.1E+1/1.2E+1);
    ctypeRNum  t360 = t170*t172*t299*(1.1E+1/1.2E+1);
    ctypeRNum  t361 = t171*t172*t299*(1.1E+1/1.2E+1);
    ctypeRNum  t362 = t137*t176*t277*(1.1E+1/2.4E+1);
    ctypeRNum  t363 = t138*t177*t277*(1.1E+1/2.4E+1);
    ctypeRNum  t364 = t139*t178*t277*(1.1E+1/2.4E+1);
    ctypeRNum  t365 = t140*t179*t279*(1.1E+1/2.4E+1);
    ctypeRNum  t366 = t141*t180*t279*(1.1E+1/2.4E+1);
    ctypeRNum  t367 = t142*t181*t279*(1.1E+1/2.4E+1);
    ctypeRNum  t368 = t143*t182*t281*(1.1E+1/2.4E+1);
    ctypeRNum  t369 = t144*t183*t281*(1.1E+1/2.4E+1);
    ctypeRNum  t370 = t145*t184*t281*(1.1E+1/2.4E+1);
    ctypeRNum  t371 = t146*t185*t283*(1.1E+1/2.4E+1);
    ctypeRNum  t372 = t147*t186*t283*(1.1E+1/2.4E+1);
    ctypeRNum  t373 = t148*t187*t283*(1.1E+1/2.4E+1);
    ctypeRNum  t374 = t149*t188*t285*(1.1E+1/2.4E+1);
    ctypeRNum  t375 = t150*t189*t285*(1.1E+1/2.4E+1);
    ctypeRNum  t376 = t151*t190*t285*(1.1E+1/2.4E+1);
    ctypeRNum  t377 = t152*t191*t287*(1.1E+1/2.4E+1);
    ctypeRNum  t378 = t153*t192*t287*(1.1E+1/2.4E+1);
    ctypeRNum  t379 = t154*t193*t287*(1.1E+1/2.4E+1);
    ctypeRNum  t380 = t155*t194*t289*(1.1E+1/2.4E+1);
    ctypeRNum  t381 = t156*t195*t289*(1.1E+1/2.4E+1);
    ctypeRNum  t382 = t157*t196*t289*(1.1E+1/2.4E+1);
    ctypeRNum  t383 = t158*t197*t291*(1.1E+1/2.4E+1);
    ctypeRNum  t384 = t159*t198*t291*(1.1E+1/2.4E+1);
    ctypeRNum  t385 = t160*t199*t291*(1.1E+1/2.4E+1);
    ctypeRNum  t386 = t161*t200*t293*(1.1E+1/2.4E+1);
    ctypeRNum  t387 = t162*t201*t293*(1.1E+1/2.4E+1);
    ctypeRNum  t388 = t163*t202*t293*(1.1E+1/2.4E+1);
    ctypeRNum  t389 = t164*t203*t295*(1.1E+1/2.4E+1);
    ctypeRNum  t390 = t165*t204*t295*(1.1E+1/2.4E+1);
    ctypeRNum  t391 = t166*t205*t295*(1.1E+1/2.4E+1);
    ctypeRNum  t392 = t167*t206*t297*(1.1E+1/2.4E+1);
    ctypeRNum  t393 = t168*t207*t297*(1.1E+1/2.4E+1);
    ctypeRNum  t394 = t169*t208*t297*(1.1E+1/2.4E+1);
    ctypeRNum  t395 = t170*t209*t299*(1.1E+1/2.4E+1);
    ctypeRNum  t396 = t171*t210*t299*(1.1E+1/2.4E+1);
    ctypeRNum  t397 = t172*t211*t299*(1.1E+1/2.4E+1);
    ctypeRNum  t398 = t259+t326;
    ctypeRNum  t399 = t260+t327;
    ctypeRNum  t400 = t261+t328;
    ctypeRNum  t401 = t314+t362+7.0E+1/3.0;
    ctypeRNum  t402 = t314+t363+7.0E+1/3.0;
    ctypeRNum  t403 = t314+t364+7.0E+1/3.0;
    ctypeRNum  t404 = t315+t365+7.0E+1/3.0;
    ctypeRNum  t405 = t315+t366+7.0E+1/3.0;
    ctypeRNum  t406 = t315+t367+7.0E+1/3.0;
    ctypeRNum  t407 = t316+t368+7.0E+1/3.0;
    ctypeRNum  t408 = t316+t369+7.0E+1/3.0;
    ctypeRNum  t409 = t316+t370+7.0E+1/3.0;
    ctypeRNum  t410 = t317+t371+7.0E+1/3.0;
    ctypeRNum  t411 = t317+t372+7.0E+1/3.0;
    ctypeRNum  t412 = t317+t373+7.0E+1/3.0;
    ctypeRNum  t413 = t318+t374+7.0E+1/3.0;
    ctypeRNum  t414 = t318+t375+7.0E+1/3.0;
    ctypeRNum  t415 = t318+t376+7.0E+1/3.0;
    ctypeRNum  t416 = t319+t377+7.0E+1/3.0;
    ctypeRNum  t417 = t319+t378+7.0E+1/3.0;
    ctypeRNum  t418 = t319+t379+7.0E+1/3.0;
    ctypeRNum  t419 = t320+t380+7.0E+1/3.0;
    ctypeRNum  t420 = t320+t381+7.0E+1/3.0;
    ctypeRNum  t421 = t320+t382+7.0E+1/3.0;
    ctypeRNum  t422 = t321+t383+7.0E+1/3.0;
    ctypeRNum  t423 = t321+t384+7.0E+1/3.0;
    ctypeRNum  t424 = t321+t385+7.0E+1/3.0;
    ctypeRNum  t425 = t322+t386+7.0E+1/3.0;
    ctypeRNum  t426 = t322+t387+7.0E+1/3.0;
    ctypeRNum  t427 = t322+t388+7.0E+1/3.0;
    ctypeRNum  t428 = t323+t389+7.0E+1/3.0;
    ctypeRNum  t429 = t323+t390+7.0E+1/3.0;
    ctypeRNum  t430 = t323+t391+7.0E+1/3.0;
    ctypeRNum  t431 = t324+t392+7.0E+1/3.0;
    ctypeRNum  t432 = t324+t393+7.0E+1/3.0;
    ctypeRNum  t433 = t324+t394+7.0E+1/3.0;
    ctypeRNum  t434 = t325+t395+7.0E+1/3.0;
    ctypeRNum  t435 = t325+t396+7.0E+1/3.0;
    ctypeRNum  t436 = t325+t397+7.0E+1/3.0;
    out[0] = t326*vec[10]+t327*vec[11]-t398*vec[4]-t399*vec[5]+t401*vec[9]-vec[3]*(t258+t314+t362+t44*t256*(1.1E+1/1.2E+1)+1.4E+2/3.0);
    out[1] = t326*vec[9]+t328*vec[11]-t398*vec[3]-t400*vec[5]+t402*vec[10]-vec[4]*(t258+t314+t363+t45*t256*(1.1E+1/1.2E+1)+1.4E+2/3.0);
    out[2] = t327*vec[9]+t328*vec[10]-t399*vec[3]-t400*vec[4]+t403*vec[11]-vec[5]*(t258+t314+t364+t46*t256*(1.1E+1/1.2E+1)+1.4E+2/3.0);
    out[3] = vec[0];
    out[4] = vec[1];
    out[5] = vec[2];
    out[6] = vec[9]*(-1.4E+2/3.0)+t302*vec[9]+t303*vec[9]+t326*vec[4]+t327*vec[5]+t329*vec[16]+t330*vec[17]+t401*vec[3]+t404*vec[15]-vec[10]*(t326+t141*t179*t279*(1.1E+1/2.4E+1))-vec[11]*(t327+t142*t179*t279*(1.1E+1/2.4E+1))-t137*t176*t277*vec[9]*(1.1E+1/2.4E+1)-t140*t179*t279*vec[9]*(1.1E+1/2.4E+1);
    out[7] = vec[10]*(-1.4E+2/3.0)+t302*vec[10]+t303*vec[10]+t326*vec[3]+t328*vec[5]+t329*vec[15]+t331*vec[17]+t402*vec[4]+t405*vec[16]-vec[9]*(t326+t140*t180*t279*(1.1E+1/2.4E+1))-vec[11]*(t328+t142*t180*t279*(1.1E+1/2.4E+1))-t138*t177*t277*vec[10]*(1.1E+1/2.4E+1)-t141*t180*t279*vec[10]*(1.1E+1/2.4E+1);
    out[8] = vec[11]*(-1.4E+2/3.0)+t302*vec[11]+t303*vec[11]+t327*vec[3]+t328*vec[4]+t330*vec[15]+t331*vec[16]+t403*vec[5]+t406*vec[17]-vec[9]*(t327+t140*t181*t279*(1.1E+1/2.4E+1))-vec[10]*(t328+t141*t181*t279*(1.1E+1/2.4E+1))-t139*t178*t277*vec[11]*(1.1E+1/2.4E+1)-t142*t181*t279*vec[11]*(1.1E+1/2.4E+1);
    out[9] = vec[6];
    out[10] = vec[7];
    out[11] = vec[8];
    out[12] = vec[15]*(-1.4E+2/3.0)+t303*vec[15]+t304*vec[15]+t329*vec[10]+t330*vec[11]+t332*vec[22]+t333*vec[23]+t404*vec[9]+t407*vec[21]-vec[16]*(t329+t144*t182*t281*(1.1E+1/2.4E+1))-vec[17]*(t330+t145*t182*t281*(1.1E+1/2.4E+1))-t140*t179*t279*vec[15]*(1.1E+1/2.4E+1)-t143*t182*t281*vec[15]*(1.1E+1/2.4E+1);
    out[13] = vec[16]*(-1.4E+2/3.0)+t303*vec[16]+t304*vec[16]+t329*vec[9]+t331*vec[11]+t332*vec[21]+t334*vec[23]+t405*vec[10]+t408*vec[22]-vec[15]*(t329+t143*t183*t281*(1.1E+1/2.4E+1))-vec[17]*(t331+t145*t183*t281*(1.1E+1/2.4E+1))-t141*t180*t279*vec[16]*(1.1E+1/2.4E+1)-t144*t183*t281*vec[16]*(1.1E+1/2.4E+1);
    out[14] = vec[17]*(-1.4E+2/3.0)+t303*vec[17]+t304*vec[17]+t330*vec[9]+t331*vec[10]+t333*vec[21]+t334*vec[22]+t406*vec[11]+t409*vec[23]-vec[15]*(t330+t143*t184*t281*(1.1E+1/2.4E+1))-vec[16]*(t331+t144*t184*t281*(1.1E+1/2.4E+1))-t142*t181*t279*vec[17]*(1.1E+1/2.4E+1)-t145*t184*t281*vec[17]*(1.1E+1/2.4E+1);
    out[15] = vec[12];
    out[16] = vec[13];
    out[17] = vec[14];
    out[18] = vec[21]*(-1.4E+2/3.0)+t304*vec[21]+t305*vec[21]+t332*vec[16]+t333*vec[17]+t335*vec[28]+t336*vec[29]+t407*vec[15]+t410*vec[27]-vec[22]*(t332+t147*t185*t283*(1.1E+1/2.4E+1))-vec[23]*(t333+t148*t185*t283*(1.1E+1/2.4E+1))-t143*t182*t281*vec[21]*(1.1E+1/2.4E+1)-t146*t185*t283*vec[21]*(1.1E+1/2.4E+1);
    out[19] = vec[22]*(-1.4E+2/3.0)+t304*vec[22]+t305*vec[22]+t332*vec[15]+t334*vec[17]+t335*vec[27]+t337*vec[29]+t408*vec[16]+t411*vec[28]-vec[21]*(t332+t146*t186*t283*(1.1E+1/2.4E+1))-vec[23]*(t334+t148*t186*t283*(1.1E+1/2.4E+1))-t144*t183*t281*vec[22]*(1.1E+1/2.4E+1)-t147*t186*t283*vec[22]*(1.1E+1/2.4E+1);
    out[20] = vec[23]*(-1.4E+2/3.0)+t304*vec[23]+t305*vec[23]+t333*vec[15]+t334*vec[16]+t336*vec[27]+t337*vec[28]+t409*vec[17]+t412*vec[29]-vec[21]*(t333+t146*t187*t283*(1.1E+1/2.4E+1))-vec[22]*(t334+t147*t187*t283*(1.1E+1/2.4E+1))-t145*t184*t281*vec[23]*(1.1E+1/2.4E+1)-t148*t187*t283*vec[23]*(1.1E+1/2.4E+1);
    out[21] = vec[18];
    out[22] = vec[19];
    out[23] = vec[20];
    out[24] = vec[27]*(-1.4E+2/3.0)+t305*vec[27]+t306*vec[27]+t335*vec[22]+t336*vec[23]+t338*vec[34]+t339*vec[35]+t410*vec[21]+t413*vec[33]-vec[28]*(t335+t150*t188*t285*(1.1E+1/2.4E+1))-vec[29]*(t336+t151*t188*t285*(1.1E+1/2.4E+1))-t146*t185*t283*vec[27]*(1.1E+1/2.4E+1)-t149*t188*t285*vec[27]*(1.1E+1/2.4E+1);
    out[25] = vec[28]*(-1.4E+2/3.0)+t305*vec[28]+t306*vec[28]+t335*vec[21]+t337*vec[23]+t338*vec[33]+t340*vec[35]+t411*vec[22]+t414*vec[34]-vec[27]*(t335+t149*t189*t285*(1.1E+1/2.4E+1))-vec[29]*(t337+t151*t189*t285*(1.1E+1/2.4E+1))-t147*t186*t283*vec[28]*(1.1E+1/2.4E+1)-t150*t189*t285*vec[28]*(1.1E+1/2.4E+1);
    out[26] = vec[29]*(-1.4E+2/3.0)+t305*vec[29]+t306*vec[29]+t336*vec[21]+t337*vec[22]+t339*vec[33]+t340*vec[34]+t412*vec[23]+t415*vec[35]-vec[27]*(t336+t149*t190*t285*(1.1E+1/2.4E+1))-vec[28]*(t337+t150*t190*t285*(1.1E+1/2.4E+1))-t148*t187*t283*vec[29]*(1.1E+1/2.4E+1)-t151*t190*t285*vec[29]*(1.1E+1/2.4E+1);
    out[27] = vec[24];
    out[28] = vec[25];
    out[29] = vec[26];
    out[30] = vec[33]*(-1.4E+2/3.0)+t306*vec[33]+t307*vec[33]+t338*vec[28]+t339*vec[29]+t341*vec[40]+t342*vec[41]+t413*vec[27]+t416*vec[39]-vec[34]*(t338+t153*t191*t287*(1.1E+1/2.4E+1))-vec[35]*(t339+t154*t191*t287*(1.1E+1/2.4E+1))-t149*t188*t285*vec[33]*(1.1E+1/2.4E+1)-t152*t191*t287*vec[33]*(1.1E+1/2.4E+1);
    out[31] = vec[34]*(-1.4E+2/3.0)+t306*vec[34]+t307*vec[34]+t338*vec[27]+t340*vec[29]+t341*vec[39]+t343*vec[41]+t414*vec[28]+t417*vec[40]-vec[33]*(t338+t152*t192*t287*(1.1E+1/2.4E+1))-vec[35]*(t340+t154*t192*t287*(1.1E+1/2.4E+1))-t150*t189*t285*vec[34]*(1.1E+1/2.4E+1)-t153*t192*t287*vec[34]*(1.1E+1/2.4E+1);
    out[32] = vec[35]*(-1.4E+2/3.0)+t306*vec[35]+t307*vec[35]+t339*vec[27]+t340*vec[28]+t342*vec[39]+t343*vec[40]+t415*vec[29]+t418*vec[41]-vec[33]*(t339+t152*t193*t287*(1.1E+1/2.4E+1))-vec[34]*(t340+t153*t193*t287*(1.1E+1/2.4E+1))-t151*t190*t285*vec[35]*(1.1E+1/2.4E+1)-t154*t193*t287*vec[35]*(1.1E+1/2.4E+1);
    out[33] = vec[30];
    out[34] = vec[31];
    out[35] = vec[32];
    out[36] = vec[39]*(-1.4E+2/3.0)+t307*vec[39]+t308*vec[39]+t341*vec[34]+t342*vec[35]+t344*vec[46]+t345*vec[47]+t416*vec[33]+t419*vec[45]-vec[40]*(t341+t156*t194*t289*(1.1E+1/2.4E+1))-vec[41]*(t342+t157*t194*t289*(1.1E+1/2.4E+1))-t152*t191*t287*vec[39]*(1.1E+1/2.4E+1)-t155*t194*t289*vec[39]*(1.1E+1/2.4E+1);
    out[37] = vec[40]*(-1.4E+2/3.0)+t307*vec[40]+t308*vec[40]+t341*vec[33]+t343*vec[35]+t344*vec[45]+t346*vec[47]+t417*vec[34]+t420*vec[46]-vec[39]*(t341+t155*t195*t289*(1.1E+1/2.4E+1))-vec[41]*(t343+t157*t195*t289*(1.1E+1/2.4E+1))-t153*t192*t287*vec[40]*(1.1E+1/2.4E+1)-t156*t195*t289*vec[40]*(1.1E+1/2.4E+1);
    out[38] = vec[41]*(-1.4E+2/3.0)+t307*vec[41]+t308*vec[41]+t342*vec[33]+t343*vec[34]+t345*vec[45]+t346*vec[46]+t418*vec[35]+t421*vec[47]-vec[39]*(t342+t155*t196*t289*(1.1E+1/2.4E+1))-vec[40]*(t343+t156*t196*t289*(1.1E+1/2.4E+1))-t154*t193*t287*vec[41]*(1.1E+1/2.4E+1)-t157*t196*t289*vec[41]*(1.1E+1/2.4E+1);
    out[39] = vec[36];
    out[40] = vec[37];
    out[41] = vec[38];
    out[42] = vec[45]*(-1.4E+2/3.0)+t308*vec[45]+t309*vec[45]+t344*vec[40]+t345*vec[41]+t347*vec[52]+t348*vec[53]+t419*vec[39]+t422*vec[51]-vec[46]*(t344+t159*t197*t291*(1.1E+1/2.4E+1))-vec[47]*(t345+t160*t197*t291*(1.1E+1/2.4E+1))-t155*t194*t289*vec[45]*(1.1E+1/2.4E+1)-t158*t197*t291*vec[45]*(1.1E+1/2.4E+1);
    out[43] = vec[46]*(-1.4E+2/3.0)+t308*vec[46]+t309*vec[46]+t344*vec[39]+t346*vec[41]+t347*vec[51]+t349*vec[53]+t420*vec[40]+t423*vec[52]-vec[45]*(t344+t158*t198*t291*(1.1E+1/2.4E+1))-vec[47]*(t346+t160*t198*t291*(1.1E+1/2.4E+1))-t156*t195*t289*vec[46]*(1.1E+1/2.4E+1)-t159*t198*t291*vec[46]*(1.1E+1/2.4E+1);
    out[44] = vec[47]*(-1.4E+2/3.0)+t308*vec[47]+t309*vec[47]+t345*vec[39]+t346*vec[40]+t348*vec[51]+t349*vec[52]+t421*vec[41]+t424*vec[53]-vec[45]*(t345+t158*t199*t291*(1.1E+1/2.4E+1))-vec[46]*(t346+t159*t199*t291*(1.1E+1/2.4E+1))-t157*t196*t289*vec[47]*(1.1E+1/2.4E+1)-t160*t199*t291*vec[47]*(1.1E+1/2.4E+1);
    out[45] = vec[42];
    out[46] = vec[43];
    out[47] = vec[44];
    out[48] = vec[51]*(-1.4E+2/3.0)+t309*vec[51]+t310*vec[51]+t347*vec[46]+t348*vec[47]+t350*vec[58]+t351*vec[59]+t422*vec[45]+t425*vec[57]-vec[52]*(t347+t162*t200*t293*(1.1E+1/2.4E+1))-vec[53]*(t348+t163*t200*t293*(1.1E+1/2.4E+1))-t158*t197*t291*vec[51]*(1.1E+1/2.4E+1)-t161*t200*t293*vec[51]*(1.1E+1/2.4E+1);
    out[49] = vec[52]*(-1.4E+2/3.0)+t309*vec[52]+t310*vec[52]+t347*vec[45]+t349*vec[47]+t350*vec[57]+t352*vec[59]+t423*vec[46]+t426*vec[58]-vec[51]*(t347+t161*t201*t293*(1.1E+1/2.4E+1))-vec[53]*(t349+t163*t201*t293*(1.1E+1/2.4E+1))-t159*t198*t291*vec[52]*(1.1E+1/2.4E+1)-t162*t201*t293*vec[52]*(1.1E+1/2.4E+1);
    out[50] = vec[53]*(-1.4E+2/3.0)+t309*vec[53]+t310*vec[53]+t348*vec[45]+t349*vec[46]+t351*vec[57]+t352*vec[58]+t424*vec[47]+t427*vec[59]-vec[51]*(t348+t161*t202*t293*(1.1E+1/2.4E+1))-vec[52]*(t349+t162*t202*t293*(1.1E+1/2.4E+1))-t160*t199*t291*vec[53]*(1.1E+1/2.4E+1)-t163*t202*t293*vec[53]*(1.1E+1/2.4E+1);
    out[51] = vec[48];
    out[52] = vec[49];
    out[53] = vec[50];
    out[54] = vec[57]*(-1.4E+2/3.0)+t310*vec[57]+t311*vec[57]+t350*vec[52]+t351*vec[53]+t353*vec[64]+t354*vec[65]+t425*vec[51]+t428*vec[63]-vec[58]*(t350+t165*t203*t295*(1.1E+1/2.4E+1))-vec[59]*(t351+t166*t203*t295*(1.1E+1/2.4E+1))-t161*t200*t293*vec[57]*(1.1E+1/2.4E+1)-t164*t203*t295*vec[57]*(1.1E+1/2.4E+1);
    out[55] = vec[58]*(-1.4E+2/3.0)+t310*vec[58]+t311*vec[58]+t350*vec[51]+t352*vec[53]+t353*vec[63]+t355*vec[65]+t426*vec[52]+t429*vec[64]-vec[57]*(t350+t164*t204*t295*(1.1E+1/2.4E+1))-vec[59]*(t352+t166*t204*t295*(1.1E+1/2.4E+1))-t162*t201*t293*vec[58]*(1.1E+1/2.4E+1)-t165*t204*t295*vec[58]*(1.1E+1/2.4E+1);
    out[56] = vec[59]*(-1.4E+2/3.0)+t310*vec[59]+t311*vec[59]+t351*vec[51]+t352*vec[52]+t354*vec[63]+t355*vec[64]+t427*vec[53]+t430*vec[65]-vec[57]*(t351+t164*t205*t295*(1.1E+1/2.4E+1))-vec[58]*(t352+t165*t205*t295*(1.1E+1/2.4E+1))-t163*t202*t293*vec[59]*(1.1E+1/2.4E+1)-t166*t205*t295*vec[59]*(1.1E+1/2.4E+1);
    out[57] = vec[54];
    out[58] = vec[55];
    out[59] = vec[56];
    out[60] = vec[63]*(-1.4E+2/3.0)+t311*vec[63]+t312*vec[63]+t353*vec[58]+t354*vec[59]+t356*vec[70]+t357*vec[71]+t428*vec[57]+t431*vec[69]-vec[64]*(t353+t168*t206*t297*(1.1E+1/2.4E+1))-vec[65]*(t354+t169*t206*t297*(1.1E+1/2.4E+1))-t164*t203*t295*vec[63]*(1.1E+1/2.4E+1)-t167*t206*t297*vec[63]*(1.1E+1/2.4E+1);
    out[61] = vec[64]*(-1.4E+2/3.0)+t311*vec[64]+t312*vec[64]+t353*vec[57]+t355*vec[59]+t356*vec[69]+t358*vec[71]+t429*vec[58]+t432*vec[70]-vec[63]*(t353+t167*t207*t297*(1.1E+1/2.4E+1))-vec[65]*(t355+t169*t207*t297*(1.1E+1/2.4E+1))-t165*t204*t295*vec[64]*(1.1E+1/2.4E+1)-t168*t207*t297*vec[64]*(1.1E+1/2.4E+1);
    out[62] = vec[65]*(-1.4E+2/3.0)+t311*vec[65]+t312*vec[65]+t354*vec[57]+t355*vec[58]+t357*vec[69]+t358*vec[70]+t430*vec[59]+t433*vec[71]-vec[63]*(t354+t167*t208*t297*(1.1E+1/2.4E+1))-vec[64]*(t355+t168*t208*t297*(1.1E+1/2.4E+1))-t166*t205*t295*vec[65]*(1.1E+1/2.4E+1)-t169*t208*t297*vec[65]*(1.1E+1/2.4E+1);
    out[63] = vec[60];
    out[64] = vec[61];
    out[65] = vec[62];
    out[66] = vec[69]*(-1.4E+2/3.0)+t312*vec[69]+t313*vec[69]+t356*vec[64]+t357*vec[65]+t359*vec[76]+t360*vec[77]+t431*vec[63]+t434*vec[75]-vec[70]*(t356+t171*t209*t299*(1.1E+1/2.4E+1))-vec[71]*(t357+t172*t209*t299*(1.1E+1/2.4E+1))-t167*t206*t297*vec[69]*(1.1E+1/2.4E+1)-t170*t209*t299*vec[69]*(1.1E+1/2.4E+1);
    out[67] = vec[70]*(-1.4E+2/3.0)+t312*vec[70]+t313*vec[70]+t356*vec[63]+t358*vec[65]+t359*vec[75]+t361*vec[77]+t432*vec[64]+t435*vec[76]-vec[69]*(t356+t170*t210*t299*(1.1E+1/2.4E+1))-vec[71]*(t358+t172*t210*t299*(1.1E+1/2.4E+1))-t168*t207*t297*vec[70]*(1.1E+1/2.4E+1)-t171*t210*t299*vec[70]*(1.1E+1/2.4E+1);
    out[68] = vec[71]*(-1.4E+2/3.0)+t312*vec[71]+t313*vec[71]+t357*vec[63]+t358*vec[64]+t360*vec[75]+t361*vec[76]+t433*vec[65]+t436*vec[77]-vec[69]*(t357+t170*t211*t299*(1.1E+1/2.4E+1))-vec[70]*(t358+t171*t211*t299*(1.1E+1/2.4E+1))-t169*t208*t297*vec[71]*(1.1E+1/2.4E+1)-t172*t211*t299*vec[71]*(1.1E+1/2.4E+1);
    out[69] = vec[66];
    out[70] = vec[67];
    out[71] = vec[68];
    out[72] = vec[75]*(-1.4E+2/3.0)+t300*vec[75]*(1.1E+1/1.2E+1)+t313*vec[75]+t359*vec[70]+t360*vec[71]+t434*vec[69]-vec[76]*(t359+t174*t212*t301*(1.1E+1/2.4E+1))-vec[77]*(t360+t175*t212*t301*(1.1E+1/2.4E+1))-t170*t209*t299*vec[75]*(1.1E+1/2.4E+1)-t173*t212*t301*vec[75]*(1.1E+1/2.4E+1);
    out[73] = vec[76]*(-1.4E+2/3.0)+t300*vec[76]*(1.1E+1/1.2E+1)+t313*vec[76]+t359*vec[69]+t361*vec[71]+t435*vec[70]-vec[75]*(t359+t173*t213*t301*(1.1E+1/2.4E+1))-vec[77]*(t361+t175*t213*t301*(1.1E+1/2.4E+1))-t171*t210*t299*vec[76]*(1.1E+1/2.4E+1)-t174*t213*t301*vec[76]*(1.1E+1/2.4E+1);
    out[74] = vec[77]*(-1.4E+2/3.0)+t300*vec[77]*(1.1E+1/1.2E+1)+t313*vec[77]+t360*vec[69]+t361*vec[70]+t436*vec[71]-vec[75]*(t360+t173*t214*t301*(1.1E+1/2.4E+1))-vec[76]*(t361+t174*t214*t301*(1.1E+1/2.4E+1))-t172*t211*t299*vec[77]*(1.1E+1/2.4E+1)-t175*t214*t301*vec[77]*(1.1E+1/2.4E+1);
    out[75] = vec[72];
    out[76] = vec[73];
    out[77] = vec[74];
    out[78] = (t301*(t47*vec[75]*1.1E+1+t275*vec[75]*2.8E+2+t129*x[72]+t130*x[72]-vec[75]*(t47+t252+t253)*1.1E+1-vec[76]*x[72]*x[79]*1.1E+1-vec[77]*x[72]*x[80]*1.1E+1))/1.2E+1-(t301*x[78]*(t129+t130+t135+t136))/1.2E+1;
    out[79] = (t301*(t48*vec[76]*1.1E+1+t275*vec[76]*2.8E+2+t128*x[73]+t130*x[73]-vec[76]*(t48+t251+t253)*1.1E+1-vec[75]*x[73]*x[78]*1.1E+1-vec[77]*x[73]*x[80]*1.1E+1))/1.2E+1-(t301*x[79]*(t128+t130+t134+t136))/1.2E+1;
    out[80] = (t301*(t49*vec[77]*1.1E+1+t275*vec[77]*2.8E+2+t128*x[74]+t129*x[74]-vec[77]*(t49+t251+t252)*1.1E+1-vec[75]*x[74]*x[78]*1.1E+1-vec[76]*x[74]*x[79]*1.1E+1))/1.2E+1-(t301*x[80]*(t128+t129+t134+t135))/1.2E+1;
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Chain_14_ProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
    out[0] = vec[78];
    out[1] = vec[79];
    out[2] = vec[80];
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void Chain_14_ProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void Chain_14_ProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*(u[0]*u[0])+pCost_[2]*(u[1]*u[1])+pCost_[2]*(u[2]*u[2])+pCost_[1]*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+pCost_[1]*(x[9]*x[9]+x[10]*x[10]+x[11]*x[11])+pCost_[1]*(x[15]*x[15]+x[16]*x[16]+x[17]*x[17])+pCost_[1]*(x[21]*x[21]+x[22]*x[22]+x[23]*x[23])+pCost_[1]*(x[27]*x[27]+x[28]*x[28]+x[29]*x[29])+pCost_[1]*(x[33]*x[33]+x[34]*x[34]+x[35]*x[35])+pCost_[1]*(x[39]*x[39]+x[40]*x[40]+x[41]*x[41])+pCost_[1]*(x[45]*x[45]+x[46]*x[46]+x[47]*x[47])+pCost_[1]*(x[51]*x[51]+x[52]*x[52]+x[53]*x[53])+pCost_[1]*(x[57]*x[57]+x[58]*x[58]+x[59]*x[59])+pCost_[1]*(x[63]*x[63]+x[64]*x[64]+x[65]*x[65])+pCost_[1]*(x[69]*x[69]+x[70]*x[70]+x[71]*x[71])+pCost_[1]*(x[75]*x[75]+x[76]*x[76]+x[77]*x[77])+pCost_[0]*(POW2(x[78]-7.0)+x[79]*x[79]+x[80]*x[80]);
}


/** Gradient dl/dx **/
void Chain_14_ProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
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
    out[66] = 0.0;
    out[67] = 0.0;
    out[68] = 0.0;
    out[69] = pCost_[1]*x[69]*2.0;
    out[70] = pCost_[1]*x[70]*2.0;
    out[71] = pCost_[1]*x[71]*2.0;
    out[72] = 0.0;
    out[73] = 0.0;
    out[74] = 0.0;
    out[75] = pCost_[1]*x[75]*2.0;
    out[76] = pCost_[1]*x[76]*2.0;
    out[77] = pCost_[1]*x[77]*2.0;
    out[78] = pCost_[0]*(x[78]-7.0)*2.0;
    out[79] = pCost_[0]*x[79]*2.0;
    out[80] = pCost_[0]*x[80]*2.0;
}


/** Gradient dl/du **/
void Chain_14_ProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    

    out[0] = pCost_[2]*u[0]*2.0;
    out[1] = pCost_[2]*u[1]*2.0;
    out[2] = pCost_[2]*u[2]*2.0;
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Chain_14_ProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Gradient dV/dx **/
void Chain_14_ProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
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
}


/** Gradient dV/dT **/
void Chain_14_ProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void Chain_14_ProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Chain_14_ProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Chain_14_ProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void Chain_14_ProblemDescription::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void Chain_14_ProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}


/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void Chain_14_ProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}