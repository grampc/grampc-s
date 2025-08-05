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
* This probfct-file describes the nonlinear chain problem with 8 chain elements from
* Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
* Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
* Process Control 22(3), 540-550 (2012)
*
*/

#include "NLChain_8_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

using namespace grampc;


Chain_8_ProblemDescription::Chain_8_ProblemDescription(const std::vector<typeRNum>& pCost)
: pCost_(pCost)
{
}


/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void Chain_8_ProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx  = 45;
    *Nu  = 3;
    *Np  = 0;
    *Nh  = 0;
    *Ng  = 0;
    *NgT = 0;
    *NhT = 0;
}


/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void Chain_8_ProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
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
    ctypeRNum  t26 = t5+x[0];
    ctypeRNum  t27 = t6+x[1];
    ctypeRNum  t28 = t7+x[2];
    ctypeRNum  t29 = t8+x[6];
    ctypeRNum  t30 = t9+x[7];
    ctypeRNum  t31 = t10+x[8];
    ctypeRNum  t32 = t11+x[12];
    ctypeRNum  t33 = t12+x[13];
    ctypeRNum  t34 = t13+x[14];
    ctypeRNum  t35 = t14+x[18];
    ctypeRNum  t36 = t15+x[19];
    ctypeRNum  t37 = t16+x[20];
    ctypeRNum  t38 = t17+x[24];
    ctypeRNum  t39 = t18+x[25];
    ctypeRNum  t40 = t19+x[26];
    ctypeRNum  t41 = t20+x[30];
    ctypeRNum  t42 = t21+x[31];
    ctypeRNum  t43 = t22+x[32];
    ctypeRNum  t44 = t23+x[36];
    ctypeRNum  t45 = t24+x[37];
    ctypeRNum  t46 = t25+x[38];
    ctypeRNum  t68 = t2+t3+t4;
    ctypeRNum  t47 = t26*t26;
    ctypeRNum  t48 = t27*t27;
    ctypeRNum  t49 = t28*t28;
    ctypeRNum  t50 = t29*t29;
    ctypeRNum  t51 = t30*t30;
    ctypeRNum  t52 = t31*t31;
    ctypeRNum  t53 = t32*t32;
    ctypeRNum  t54 = t33*t33;
    ctypeRNum  t55 = t34*t34;
    ctypeRNum  t56 = t35*t35;
    ctypeRNum  t57 = t36*t36;
    ctypeRNum  t58 = t37*t37;
    ctypeRNum  t59 = t38*t38;
    ctypeRNum  t60 = t39*t39;
    ctypeRNum  t61 = t40*t40;
    ctypeRNum  t62 = t41*t41;
    ctypeRNum  t63 = t42*t42;
    ctypeRNum  t64 = t43*t43;
    ctypeRNum  t65 = t44*t44;
    ctypeRNum  t66 = t45*t45;
    ctypeRNum  t67 = t46*t46;
    ctypeRNum  t69 = 1.0/sqrt(t68);
    ctypeRNum  t70 = t69*(1.1E+1/4.0E+2);
    ctypeRNum  t71 = t47+t48+t49;
    ctypeRNum  t72 = t50+t51+t52;
    ctypeRNum  t73 = t53+t54+t55;
    ctypeRNum  t74 = t56+t57+t58;
    ctypeRNum  t75 = t59+t60+t61;
    ctypeRNum  t76 = t62+t63+t64;
    ctypeRNum  t77 = t65+t66+t67;
    ctypeRNum  t78 = t70-2.0/5.0;
    ctypeRNum  t79 = 1.0/sqrt(t71);
    ctypeRNum  t80 = 1.0/sqrt(t72);
    ctypeRNum  t81 = 1.0/sqrt(t73);
    ctypeRNum  t82 = 1.0/sqrt(t74);
    ctypeRNum  t83 = 1.0/sqrt(t75);
    ctypeRNum  t84 = 1.0/sqrt(t76);
    ctypeRNum  t85 = 1.0/sqrt(t77);
    ctypeRNum  t86 = t79*(1.1E+1/4.0E+2);
    ctypeRNum  t87 = t80*(1.1E+1/4.0E+2);
    ctypeRNum  t88 = t81*(1.1E+1/4.0E+2);
    ctypeRNum  t89 = t82*(1.1E+1/4.0E+2);
    ctypeRNum  t90 = t83*(1.1E+1/4.0E+2);
    ctypeRNum  t91 = t84*(1.1E+1/4.0E+2);
    ctypeRNum  t92 = t85*(1.1E+1/4.0E+2);
    ctypeRNum  t93 = t86-2.0/5.0;
    ctypeRNum  t94 = t87-2.0/5.0;
    ctypeRNum  t95 = t88-2.0/5.0;
    ctypeRNum  t96 = t89-2.0/5.0;
    ctypeRNum  t97 = t90-2.0/5.0;
    ctypeRNum  t98 = t91-2.0/5.0;
    ctypeRNum  t99 = t92-2.0/5.0;
    ctypeRNum  t100 = t26*t93*(1.0E+2/3.0);
    ctypeRNum  t101 = t27*t93*(1.0E+2/3.0);
    ctypeRNum  t102 = t28*t93*(1.0E+2/3.0);
    ctypeRNum  t103 = t29*t94*(1.0E+2/3.0);
    ctypeRNum  t104 = t30*t94*(1.0E+2/3.0);
    ctypeRNum  t105 = t31*t94*(1.0E+2/3.0);
    ctypeRNum  t106 = t32*t95*(1.0E+2/3.0);
    ctypeRNum  t107 = t33*t95*(1.0E+2/3.0);
    ctypeRNum  t108 = t34*t95*(1.0E+2/3.0);
    ctypeRNum  t109 = t35*t96*(1.0E+2/3.0);
    ctypeRNum  t110 = t36*t96*(1.0E+2/3.0);
    ctypeRNum  t111 = t37*t96*(1.0E+2/3.0);
    ctypeRNum  t112 = t38*t97*(1.0E+2/3.0);
    ctypeRNum  t113 = t39*t97*(1.0E+2/3.0);
    ctypeRNum  t114 = t40*t97*(1.0E+2/3.0);
    ctypeRNum  t115 = t41*t98*(1.0E+2/3.0);
    ctypeRNum  t116 = t42*t98*(1.0E+2/3.0);
    ctypeRNum  t117 = t43*t98*(1.0E+2/3.0);
    out[0] = x[3];
    out[1] = x[4];
    out[2] = x[5];
    out[3] = t100+t78*x[0]*(1.0E+2/3.0);
    out[4] = t101+t78*x[1]*(1.0E+2/3.0);
    out[5] = t102+t78*x[2]*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[6] = x[9];
    out[7] = x[10];
    out[8] = x[11];
    out[9] = -t100+t103;
    out[10] = -t101+t104;
    out[11] = -t102+t105-9.81E+2/1.0E+2;
    out[12] = x[15];
    out[13] = x[16];
    out[14] = x[17];
    out[15] = -t103+t106;
    out[16] = -t104+t107;
    out[17] = -t105+t108-9.81E+2/1.0E+2;
    out[18] = x[21];
    out[19] = x[22];
    out[20] = x[23];
    out[21] = -t106+t109;
    out[22] = -t107+t110;
    out[23] = -t108+t111-9.81E+2/1.0E+2;
    out[24] = x[27];
    out[25] = x[28];
    out[26] = x[29];
    out[27] = -t109+t112;
    out[28] = -t110+t113;
    out[29] = -t111+t114-9.81E+2/1.0E+2;
    out[30] = x[33];
    out[31] = x[34];
    out[32] = x[35];
    out[33] = -t112+t115;
    out[34] = -t113+t116;
    out[35] = -t114+t117-9.81E+2/1.0E+2;
    out[36] = x[39];
    out[37] = x[40];
    out[38] = x[41];
    out[39] = -t115+t44*t99*(1.0E+2/3.0);
    out[40] = -t116+t45*t99*(1.0E+2/3.0);
    out[41] = -t117+t46*t99*(1.0E+2/3.0)-9.81E+2/1.0E+2;
    out[42] = u[0];
    out[43] = u[1];
    out[44] = u[2];
}


/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Chain_8_ProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
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
    ctypeRNum  t26 = x[0]*x[0];
    ctypeRNum  t27 = x[1]*x[1];
    ctypeRNum  t28 = x[2]*x[2];
    ctypeRNum  t29 = x[36]*x[36];
    ctypeRNum  t30 = x[37]*x[37];
    ctypeRNum  t31 = x[38]*x[38];
    ctypeRNum  t32 = -x[6];
    ctypeRNum  t34 = -x[7];
    ctypeRNum  t36 = -x[8];
    ctypeRNum  t38 = -x[12];
    ctypeRNum  t40 = -x[13];
    ctypeRNum  t42 = -x[14];
    ctypeRNum  t44 = -x[18];
    ctypeRNum  t46 = -x[19];
    ctypeRNum  t48 = -x[20];
    ctypeRNum  t50 = -x[24];
    ctypeRNum  t52 = -x[25];
    ctypeRNum  t54 = -x[26];
    ctypeRNum  t56 = -x[30];
    ctypeRNum  t58 = -x[31];
    ctypeRNum  t60 = -x[32];
    ctypeRNum  t62 = -x[36];
    ctypeRNum  t64 = -x[37];
    ctypeRNum  t66 = -x[38];
    ctypeRNum  t68 = -x[42];
    ctypeRNum  t70 = -x[43];
    ctypeRNum  t72 = -x[44];
    ctypeRNum  t74 = vec[39]*x[36]*1.1E+1;
    ctypeRNum  t75 = vec[40]*x[37]*1.1E+1;
    ctypeRNum  t76 = vec[41]*x[38]*1.1E+1;
    ctypeRNum  t77 = vec[39]*x[42]*1.1E+1;
    ctypeRNum  t78 = vec[40]*x[43]*1.1E+1;
    ctypeRNum  t79 = vec[41]*x[44]*1.1E+1;
    ctypeRNum  t33 = -t5;
    ctypeRNum  t35 = -t6;
    ctypeRNum  t37 = -t7;
    ctypeRNum  t39 = -t8;
    ctypeRNum  t41 = -t9;
    ctypeRNum  t43 = -t10;
    ctypeRNum  t45 = -t11;
    ctypeRNum  t47 = -t12;
    ctypeRNum  t49 = -t13;
    ctypeRNum  t51 = -t14;
    ctypeRNum  t53 = -t15;
    ctypeRNum  t55 = -t16;
    ctypeRNum  t57 = -t17;
    ctypeRNum  t59 = -t18;
    ctypeRNum  t61 = -t19;
    ctypeRNum  t63 = -t20;
    ctypeRNum  t65 = -t21;
    ctypeRNum  t67 = -t22;
    ctypeRNum  t69 = -t23;
    ctypeRNum  t71 = -t24;
    ctypeRNum  t73 = -t25;
    ctypeRNum  t80 = -t77;
    ctypeRNum  t81 = -t78;
    ctypeRNum  t82 = -t79;
    ctypeRNum  t83 = t32+x[0];
    ctypeRNum  t84 = t34+x[1];
    ctypeRNum  t85 = t36+x[2];
    ctypeRNum  t86 = t38+x[6];
    ctypeRNum  t87 = t40+x[7];
    ctypeRNum  t88 = t42+x[8];
    ctypeRNum  t89 = t44+x[12];
    ctypeRNum  t90 = t46+x[13];
    ctypeRNum  t91 = t48+x[14];
    ctypeRNum  t92 = t50+x[18];
    ctypeRNum  t93 = t52+x[19];
    ctypeRNum  t94 = t54+x[20];
    ctypeRNum  t95 = t56+x[24];
    ctypeRNum  t96 = t58+x[25];
    ctypeRNum  t97 = t60+x[26];
    ctypeRNum  t98 = t62+x[30];
    ctypeRNum  t99 = t64+x[31];
    ctypeRNum  t100 = t66+x[32];
    ctypeRNum  t101 = t68+x[36];
    ctypeRNum  t102 = t70+x[37];
    ctypeRNum  t103 = t72+x[38];
    ctypeRNum  t146 = t26+t27+t28;
    ctypeRNum  t104 = t2+t33;
    ctypeRNum  t105 = t3+t35;
    ctypeRNum  t106 = t4+t37;
    ctypeRNum  t107 = t5+t39;
    ctypeRNum  t108 = t6+t41;
    ctypeRNum  t109 = t7+t43;
    ctypeRNum  t110 = t8+t45;
    ctypeRNum  t111 = t9+t47;
    ctypeRNum  t112 = t10+t49;
    ctypeRNum  t113 = t11+t51;
    ctypeRNum  t114 = t12+t53;
    ctypeRNum  t115 = t13+t55;
    ctypeRNum  t116 = t14+t57;
    ctypeRNum  t117 = t15+t59;
    ctypeRNum  t118 = t16+t61;
    ctypeRNum  t119 = t17+t63;
    ctypeRNum  t120 = t18+t65;
    ctypeRNum  t121 = t19+t67;
    ctypeRNum  t122 = t20+t69;
    ctypeRNum  t123 = t21+t71;
    ctypeRNum  t124 = t22+t73;
    ctypeRNum  t125 = t83*t83;
    ctypeRNum  t126 = t84*t84;
    ctypeRNum  t127 = t85*t85;
    ctypeRNum  t128 = t86*t86;
    ctypeRNum  t129 = t87*t87;
    ctypeRNum  t130 = t88*t88;
    ctypeRNum  t131 = t89*t89;
    ctypeRNum  t132 = t90*t90;
    ctypeRNum  t133 = t91*t91;
    ctypeRNum  t134 = t92*t92;
    ctypeRNum  t135 = t93*t93;
    ctypeRNum  t136 = t94*t94;
    ctypeRNum  t137 = t95*t95;
    ctypeRNum  t138 = t96*t96;
    ctypeRNum  t139 = t97*t97;
    ctypeRNum  t140 = t98*t98;
    ctypeRNum  t141 = t99*t99;
    ctypeRNum  t142 = t100*t100;
    ctypeRNum  t143 = t101*t101;
    ctypeRNum  t144 = t102*t102;
    ctypeRNum  t145 = t103*t103;
    ctypeRNum  t147 = 1.0/sqrt(t146);
    ctypeRNum  t148 = t147*t147*t147;
    ctypeRNum  t149 = t147*(1.1E+1/1.2E+1);
    ctypeRNum  t154 = t125+t126+t127;
    ctypeRNum  t155 = t128+t129+t130;
    ctypeRNum  t156 = t131+t132+t133;
    ctypeRNum  t157 = t134+t135+t136;
    ctypeRNum  t158 = t137+t138+t139;
    ctypeRNum  t159 = t140+t141+t142;
    ctypeRNum  t160 = t143+t144+t145;
    ctypeRNum  t150 = -t149;
    ctypeRNum  t151 = t148*x[0]*x[1]*(1.1E+1/1.2E+1);
    ctypeRNum  t152 = t148*x[0]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t153 = t148*x[1]*x[2]*(1.1E+1/1.2E+1);
    ctypeRNum  t161 = pow(t160,3.0/2.0);
    ctypeRNum  t162 = 1.0/sqrt(t154);
    ctypeRNum  t164 = 1.0/sqrt(t155);
    ctypeRNum  t166 = 1.0/sqrt(t156);
    ctypeRNum  t168 = 1.0/sqrt(t157);
    ctypeRNum  t170 = 1.0/sqrt(t158);
    ctypeRNum  t172 = 1.0/sqrt(t159);
    ctypeRNum  t174 = 1.0/sqrt(t160);
    ctypeRNum  t163 = t162*t162*t162;
    ctypeRNum  t165 = t164*t164*t164;
    ctypeRNum  t167 = t166*t166*t166;
    ctypeRNum  t169 = t168*t168*t168;
    ctypeRNum  t171 = t170*t170*t170;
    ctypeRNum  t173 = t172*t172*t172;
    ctypeRNum  t175 = 1.0/t161;
    ctypeRNum  t176 = t162*(1.1E+1/1.2E+1);
    ctypeRNum  t177 = t164*(1.1E+1/1.2E+1);
    ctypeRNum  t178 = t166*(1.1E+1/1.2E+1);
    ctypeRNum  t179 = t168*(1.1E+1/1.2E+1);
    ctypeRNum  t180 = t170*(1.1E+1/1.2E+1);
    ctypeRNum  t181 = t172*(1.1E+1/1.2E+1);
    ctypeRNum  t182 = -t176;
    ctypeRNum  t183 = -t177;
    ctypeRNum  t184 = -t178;
    ctypeRNum  t185 = -t179;
    ctypeRNum  t186 = -t180;
    ctypeRNum  t187 = -t181;
    ctypeRNum  t188 = t83*t84*t163*(1.1E+1/1.2E+1);
    ctypeRNum  t189 = t83*t85*t163*(1.1E+1/1.2E+1);
    ctypeRNum  t190 = t84*t85*t163*(1.1E+1/1.2E+1);
    ctypeRNum  t191 = t86*t87*t165*(1.1E+1/1.2E+1);
    ctypeRNum  t192 = t86*t88*t165*(1.1E+1/1.2E+1);
    ctypeRNum  t193 = t87*t88*t165*(1.1E+1/1.2E+1);
    ctypeRNum  t194 = t89*t90*t167*(1.1E+1/1.2E+1);
    ctypeRNum  t195 = t89*t91*t167*(1.1E+1/1.2E+1);
    ctypeRNum  t196 = t90*t91*t167*(1.1E+1/1.2E+1);
    ctypeRNum  t197 = t92*t93*t169*(1.1E+1/1.2E+1);
    ctypeRNum  t198 = t92*t94*t169*(1.1E+1/1.2E+1);
    ctypeRNum  t199 = t93*t94*t169*(1.1E+1/1.2E+1);
    ctypeRNum  t200 = t95*t96*t171*(1.1E+1/1.2E+1);
    ctypeRNum  t201 = t95*t97*t171*(1.1E+1/1.2E+1);
    ctypeRNum  t202 = t96*t97*t171*(1.1E+1/1.2E+1);
    ctypeRNum  t203 = t98*t99*t173*(1.1E+1/1.2E+1);
    ctypeRNum  t204 = t98*t100*t173*(1.1E+1/1.2E+1);
    ctypeRNum  t205 = t99*t100*t173*(1.1E+1/1.2E+1);
    ctypeRNum  t206 = t83*t104*t163*(1.1E+1/2.4E+1);
    ctypeRNum  t207 = t84*t105*t163*(1.1E+1/2.4E+1);
    ctypeRNum  t208 = t85*t106*t163*(1.1E+1/2.4E+1);
    ctypeRNum  t209 = t86*t107*t165*(1.1E+1/2.4E+1);
    ctypeRNum  t210 = t87*t108*t165*(1.1E+1/2.4E+1);
    ctypeRNum  t211 = t88*t109*t165*(1.1E+1/2.4E+1);
    ctypeRNum  t212 = t89*t110*t167*(1.1E+1/2.4E+1);
    ctypeRNum  t213 = t90*t111*t167*(1.1E+1/2.4E+1);
    ctypeRNum  t214 = t91*t112*t167*(1.1E+1/2.4E+1);
    ctypeRNum  t215 = t92*t113*t169*(1.1E+1/2.4E+1);
    ctypeRNum  t216 = t93*t114*t169*(1.1E+1/2.4E+1);
    ctypeRNum  t217 = t94*t115*t169*(1.1E+1/2.4E+1);
    ctypeRNum  t218 = t95*t116*t171*(1.1E+1/2.4E+1);
    ctypeRNum  t219 = t96*t117*t171*(1.1E+1/2.4E+1);
    ctypeRNum  t220 = t97*t118*t171*(1.1E+1/2.4E+1);
    ctypeRNum  t221 = t98*t119*t173*(1.1E+1/2.4E+1);
    ctypeRNum  t222 = t99*t120*t173*(1.1E+1/2.4E+1);
    ctypeRNum  t223 = t100*t121*t173*(1.1E+1/2.4E+1);
    ctypeRNum  t224 = t151+t188;
    ctypeRNum  t225 = t152+t189;
    ctypeRNum  t226 = t153+t190;
    ctypeRNum  t227 = t182+t206+4.0E+1/3.0;
    ctypeRNum  t228 = t182+t207+4.0E+1/3.0;
    ctypeRNum  t229 = t182+t208+4.0E+1/3.0;
    ctypeRNum  t230 = t183+t209+4.0E+1/3.0;
    ctypeRNum  t231 = t183+t210+4.0E+1/3.0;
    ctypeRNum  t232 = t183+t211+4.0E+1/3.0;
    ctypeRNum  t233 = t184+t212+4.0E+1/3.0;
    ctypeRNum  t234 = t184+t213+4.0E+1/3.0;
    ctypeRNum  t235 = t184+t214+4.0E+1/3.0;
    ctypeRNum  t236 = t185+t215+4.0E+1/3.0;
    ctypeRNum  t237 = t185+t216+4.0E+1/3.0;
    ctypeRNum  t238 = t185+t217+4.0E+1/3.0;
    ctypeRNum  t239 = t186+t218+4.0E+1/3.0;
    ctypeRNum  t240 = t186+t219+4.0E+1/3.0;
    ctypeRNum  t241 = t186+t220+4.0E+1/3.0;
    ctypeRNum  t242 = t187+t221+4.0E+1/3.0;
    ctypeRNum  t243 = t187+t222+4.0E+1/3.0;
    ctypeRNum  t244 = t187+t223+4.0E+1/3.0;
    out[0] = t188*vec[10]+t189*vec[11]-t224*vec[4]-t225*vec[5]+t227*vec[9]-vec[3]*(t150+t182+t206+t26*t148*(1.1E+1/1.2E+1)+8.0E+1/3.0);
    out[1] = t188*vec[9]+t190*vec[11]-t224*vec[3]-t226*vec[5]+t228*vec[10]-vec[4]*(t150+t182+t207+t27*t148*(1.1E+1/1.2E+1)+8.0E+1/3.0);
    out[2] = t189*vec[9]+t190*vec[10]-t225*vec[3]-t226*vec[4]+t229*vec[11]-vec[5]*(t150+t182+t208+t28*t148*(1.1E+1/1.2E+1)+8.0E+1/3.0);
    out[3] = vec[0];
    out[4] = vec[1];
    out[5] = vec[2];
    out[6] = vec[9]*(-8.0E+1/3.0)+t176*vec[9]+t177*vec[9]+t188*vec[4]+t189*vec[5]+t191*vec[16]+t192*vec[17]+t227*vec[3]+t230*vec[15]-vec[10]*(t188+t87*t107*t165*(1.1E+1/2.4E+1))-vec[11]*(t189+t88*t107*t165*(1.1E+1/2.4E+1))-t83*t104*t163*vec[9]*(1.1E+1/2.4E+1)-t86*t107*t165*vec[9]*(1.1E+1/2.4E+1);
    out[7] = vec[10]*(-8.0E+1/3.0)+t176*vec[10]+t177*vec[10]+t188*vec[3]+t190*vec[5]+t191*vec[15]+t193*vec[17]+t228*vec[4]+t231*vec[16]-vec[9]*(t188+t86*t108*t165*(1.1E+1/2.4E+1))-vec[11]*(t190+t88*t108*t165*(1.1E+1/2.4E+1))-t84*t105*t163*vec[10]*(1.1E+1/2.4E+1)-t87*t108*t165*vec[10]*(1.1E+1/2.4E+1);
    out[8] = vec[11]*(-8.0E+1/3.0)+t176*vec[11]+t177*vec[11]+t189*vec[3]+t190*vec[4]+t192*vec[15]+t193*vec[16]+t229*vec[5]+t232*vec[17]-vec[9]*(t189+t86*t109*t165*(1.1E+1/2.4E+1))-vec[10]*(t190+t87*t109*t165*(1.1E+1/2.4E+1))-t85*t106*t163*vec[11]*(1.1E+1/2.4E+1)-t88*t109*t165*vec[11]*(1.1E+1/2.4E+1);
    out[9] = vec[6];
    out[10] = vec[7];
    out[11] = vec[8];
    out[12] = vec[15]*(-8.0E+1/3.0)+t177*vec[15]+t178*vec[15]+t191*vec[10]+t192*vec[11]+t194*vec[22]+t195*vec[23]+t230*vec[9]+t233*vec[21]-vec[16]*(t191+t90*t110*t167*(1.1E+1/2.4E+1))-vec[17]*(t192+t91*t110*t167*(1.1E+1/2.4E+1))-t86*t107*t165*vec[15]*(1.1E+1/2.4E+1)-t89*t110*t167*vec[15]*(1.1E+1/2.4E+1);
    out[13] = vec[16]*(-8.0E+1/3.0)+t177*vec[16]+t178*vec[16]+t191*vec[9]+t193*vec[11]+t194*vec[21]+t196*vec[23]+t231*vec[10]+t234*vec[22]-vec[15]*(t191+t89*t111*t167*(1.1E+1/2.4E+1))-vec[17]*(t193+t91*t111*t167*(1.1E+1/2.4E+1))-t87*t108*t165*vec[16]*(1.1E+1/2.4E+1)-t90*t111*t167*vec[16]*(1.1E+1/2.4E+1);
    out[14] = vec[17]*(-8.0E+1/3.0)+t177*vec[17]+t178*vec[17]+t192*vec[9]+t193*vec[10]+t195*vec[21]+t196*vec[22]+t232*vec[11]+t235*vec[23]-vec[15]*(t192+t89*t112*t167*(1.1E+1/2.4E+1))-vec[16]*(t193+t90*t112*t167*(1.1E+1/2.4E+1))-t88*t109*t165*vec[17]*(1.1E+1/2.4E+1)-t91*t112*t167*vec[17]*(1.1E+1/2.4E+1);
    out[15] = vec[12];
    out[16] = vec[13];
    out[17] = vec[14];
    out[18] = vec[21]*(-8.0E+1/3.0)+t178*vec[21]+t179*vec[21]+t194*vec[16]+t195*vec[17]+t197*vec[28]+t198*vec[29]+t233*vec[15]+t236*vec[27]-vec[22]*(t194+t93*t113*t169*(1.1E+1/2.4E+1))-vec[23]*(t195+t94*t113*t169*(1.1E+1/2.4E+1))-t89*t110*t167*vec[21]*(1.1E+1/2.4E+1)-t92*t113*t169*vec[21]*(1.1E+1/2.4E+1);
    out[19] = vec[22]*(-8.0E+1/3.0)+t178*vec[22]+t179*vec[22]+t194*vec[15]+t196*vec[17]+t197*vec[27]+t199*vec[29]+t234*vec[16]+t237*vec[28]-vec[21]*(t194+t92*t114*t169*(1.1E+1/2.4E+1))-vec[23]*(t196+t94*t114*t169*(1.1E+1/2.4E+1))-t90*t111*t167*vec[22]*(1.1E+1/2.4E+1)-t93*t114*t169*vec[22]*(1.1E+1/2.4E+1);
    out[20] = vec[23]*(-8.0E+1/3.0)+t178*vec[23]+t179*vec[23]+t195*vec[15]+t196*vec[16]+t198*vec[27]+t199*vec[28]+t235*vec[17]+t238*vec[29]-vec[21]*(t195+t92*t115*t169*(1.1E+1/2.4E+1))-vec[22]*(t196+t93*t115*t169*(1.1E+1/2.4E+1))-t91*t112*t167*vec[23]*(1.1E+1/2.4E+1)-t94*t115*t169*vec[23]*(1.1E+1/2.4E+1);
    out[21] = vec[18];
    out[22] = vec[19];
    out[23] = vec[20];
    out[24] = vec[27]*(-8.0E+1/3.0)+t179*vec[27]+t180*vec[27]+t197*vec[22]+t198*vec[23]+t200*vec[34]+t201*vec[35]+t236*vec[21]+t239*vec[33]-vec[28]*(t197+t96*t116*t171*(1.1E+1/2.4E+1))-vec[29]*(t198+t97*t116*t171*(1.1E+1/2.4E+1))-t92*t113*t169*vec[27]*(1.1E+1/2.4E+1)-t95*t116*t171*vec[27]*(1.1E+1/2.4E+1);
    out[25] = vec[28]*(-8.0E+1/3.0)+t179*vec[28]+t180*vec[28]+t197*vec[21]+t199*vec[23]+t200*vec[33]+t202*vec[35]+t237*vec[22]+t240*vec[34]-vec[27]*(t197+t95*t117*t171*(1.1E+1/2.4E+1))-vec[29]*(t199+t97*t117*t171*(1.1E+1/2.4E+1))-t93*t114*t169*vec[28]*(1.1E+1/2.4E+1)-t96*t117*t171*vec[28]*(1.1E+1/2.4E+1);
    out[26] = vec[29]*(-8.0E+1/3.0)+t179*vec[29]+t180*vec[29]+t198*vec[21]+t199*vec[22]+t201*vec[33]+t202*vec[34]+t238*vec[23]+t241*vec[35]-vec[27]*(t198+t95*t118*t171*(1.1E+1/2.4E+1))-vec[28]*(t199+t96*t118*t171*(1.1E+1/2.4E+1))-t94*t115*t169*vec[29]*(1.1E+1/2.4E+1)-t97*t118*t171*vec[29]*(1.1E+1/2.4E+1);
    out[27] = vec[24];
    out[28] = vec[25];
    out[29] = vec[26];
    out[30] = vec[33]*(-8.0E+1/3.0)+t180*vec[33]+t181*vec[33]+t200*vec[28]+t201*vec[29]+t203*vec[40]+t204*vec[41]+t239*vec[27]+t242*vec[39]-vec[34]*(t200+t99*t119*t173*(1.1E+1/2.4E+1))-vec[35]*(t201+t100*t119*t173*(1.1E+1/2.4E+1))-t95*t116*t171*vec[33]*(1.1E+1/2.4E+1)-t98*t119*t173*vec[33]*(1.1E+1/2.4E+1);
    out[31] = vec[34]*(-8.0E+1/3.0)+t180*vec[34]+t181*vec[34]+t200*vec[27]+t202*vec[29]+t203*vec[39]+t205*vec[41]+t240*vec[28]+t243*vec[40]-vec[33]*(t200+t98*t120*t173*(1.1E+1/2.4E+1))-vec[35]*(t202+t100*t120*t173*(1.1E+1/2.4E+1))-t96*t117*t171*vec[34]*(1.1E+1/2.4E+1)-t99*t120*t173*vec[34]*(1.1E+1/2.4E+1);
    out[32] = vec[35]*(-8.0E+1/3.0)+t180*vec[35]+t181*vec[35]+t201*vec[27]+t202*vec[28]+t204*vec[39]+t205*vec[40]+t241*vec[29]+t244*vec[41]-vec[33]*(t201+t98*t121*t173*(1.1E+1/2.4E+1))-vec[34]*(t202+t99*t121*t173*(1.1E+1/2.4E+1))-t97*t118*t171*vec[35]*(1.1E+1/2.4E+1)-t100*t121*t173*vec[35]*(1.1E+1/2.4E+1);
    out[33] = vec[30];
    out[34] = vec[31];
    out[35] = vec[32];
    out[36] = vec[39]*(-8.0E+1/3.0)+t174*vec[39]*(1.1E+1/1.2E+1)+t181*vec[39]+t203*vec[34]+t204*vec[35]+t242*vec[33]-vec[40]*(t203+t102*t122*t175*(1.1E+1/2.4E+1))-vec[41]*(t204+t103*t122*t175*(1.1E+1/2.4E+1))-t98*t119*t173*vec[39]*(1.1E+1/2.4E+1)-t101*t122*t175*vec[39]*(1.1E+1/2.4E+1);
    out[37] = vec[40]*(-8.0E+1/3.0)+t174*vec[40]*(1.1E+1/1.2E+1)+t181*vec[40]+t203*vec[33]+t205*vec[35]+t243*vec[34]-vec[39]*(t203+t101*t123*t175*(1.1E+1/2.4E+1))-vec[41]*(t205+t103*t123*t175*(1.1E+1/2.4E+1))-t99*t120*t173*vec[40]*(1.1E+1/2.4E+1)-t102*t123*t175*vec[40]*(1.1E+1/2.4E+1);
    out[38] = vec[41]*(-8.0E+1/3.0)+t174*vec[41]*(1.1E+1/1.2E+1)+t181*vec[41]+t204*vec[33]+t205*vec[34]+t244*vec[35]-vec[39]*(t204+t101*t124*t175*(1.1E+1/2.4E+1))-vec[40]*(t205+t102*t124*t175*(1.1E+1/2.4E+1))-t100*t121*t173*vec[41]*(1.1E+1/2.4E+1)-t103*t124*t175*vec[41]*(1.1E+1/2.4E+1);
    out[39] = vec[36];
    out[40] = vec[37];
    out[41] = vec[38];
    out[42] = (t175*(t29*vec[39]*1.1E+1+t161*vec[39]*1.6E+2+t75*x[36]+t76*x[36]-vec[39]*(t29+t144+t145)*1.1E+1-vec[40]*x[36]*x[43]*1.1E+1-vec[41]*x[36]*x[44]*1.1E+1))/1.2E+1-(t175*x[42]*(t75+t76+t81+t82))/1.2E+1;
    out[43] = (t175*(t30*vec[40]*1.1E+1+t161*vec[40]*1.6E+2+t74*x[37]+t76*x[37]-vec[40]*(t30+t143+t145)*1.1E+1-vec[39]*x[37]*x[42]*1.1E+1-vec[41]*x[37]*x[44]*1.1E+1))/1.2E+1-(t175*x[43]*(t74+t76+t80+t82))/1.2E+1;
    out[44] = (t175*(t31*vec[41]*1.1E+1+t161*vec[41]*1.6E+2+t74*x[38]+t75*x[38]-vec[41]*(t31+t143+t144)*1.1E+1-vec[39]*x[38]*x[42]*1.1E+1-vec[40]*x[38]*x[43]*1.1E+1))/1.2E+1-(t175*x[44]*(t74+t75+t80+t81))/1.2E+1;
}


/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Chain_8_ProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
    out[0] = vec[42];
    out[1] = vec[43];
    out[2] = vec[44];
}


/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void Chain_8_ProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void Chain_8_ProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    

    out[0] = pCost_[2]*(u[0]*u[0])+pCost_[2]*(u[1]*u[1])+pCost_[2]*(u[2]*u[2])+pCost_[1]*(x[3]*x[3]+x[4]*x[4]+x[5]*x[5])+pCost_[1]*(x[9]*x[9]+x[10]*x[10]+x[11]*x[11])+pCost_[1]*(x[15]*x[15]+x[16]*x[16]+x[17]*x[17])+pCost_[1]*(x[21]*x[21]+x[22]*x[22]+x[23]*x[23])+pCost_[1]*(x[27]*x[27]+x[28]*x[28]+x[29]*x[29])+pCost_[1]*(x[33]*x[33]+x[34]*x[34]+x[35]*x[35])+pCost_[1]*(x[39]*x[39]+x[40]*x[40]+x[41]*x[41])+pCost_[0]*(POW2(x[42]-4.0)+x[43]*x[43]+x[44]*x[44]);
}


/** Gradient dl/dx **/
void Chain_8_ProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
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
    out[42] = pCost_[0]*(x[42]-4.0)*2.0;
    out[43] = pCost_[0]*x[43]*2.0;
    out[44] = pCost_[0]*x[44]*2.0;
}


/** Gradient dl/du **/
void Chain_8_ProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes)
{
    

    out[0] = pCost_[2]*u[0]*2.0;
    out[1] = pCost_[2]*u[1]*2.0;
    out[2] = pCost_[2]*u[2]*2.0;
}



/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Chain_8_ProblemDescription::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 0.0;
}


/** Gradient dV/dx **/
void Chain_8_ProblemDescription::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
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
}


/** Gradient dV/dT **/
void Chain_8_ProblemDescription::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes)
{
    out[0] = 0.0;
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void Chain_8_ProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
}


/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Chain_8_ProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
}


/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Chain_8_ProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec)
{
}



/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void Chain_8_ProblemDescription::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p)
{
}


/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void Chain_8_ProblemDescription::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
{
}


/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void Chain_8_ProblemDescription::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec)
{
}