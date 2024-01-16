#include "vehicle_problem_description.hpp"

VehicleProblemDescription::VehicleProblemDescription(const std::vector<typeRNum> &pSys,
                                                     const std::vector<typeRNum> &pCost,
                                                     const std::vector<typeRNum> &pCon)
    : pSys_(pSys), pCost_(pCost), pCon_(pCon)
{
    F_f_z_ = pSys_[0] * pSys_[1] * pSys_[4] / (pSys_[3] + pSys_[4]);
    F_r_z_ = pSys_[0] * pSys_[1] * pSys_[3] / (pSys_[3] + pSys_[4]);
}

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void VehicleProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 6;
    *Nu = 1;
    *Np = 4;
    *Nh = 3;
    *Ng = 0;
    *NgT = 0;
    *NhT = 0;
}

/** System function f(t,x,u,p,userparam)
    ------------------------------------ **/
void VehicleProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    ctypeRNum t2 = cos(x[2]);
    ctypeRNum t3 = cos(x[5]);
    ctypeRNum t4 = sin(x[2]);
    ctypeRNum t5 = pSys_[3] + pSys_[4];
    ctypeRNum t6 = pSys_[3] * x[1];
    ctypeRNum t7 = pSys_[4] * x[1];
    ctypeRNum t8 = x[0] + x[2];
    ctypeRNum t9 = 1.0 / 3.141592653589793;
    ctypeRNum t11 = 1.0 / pSys_[7];
    ctypeRNum t10 = pSys_[7] * t4;
    ctypeRNum t12 = 1.0 / t2;
    ctypeRNum t13 = 1.0 / t5;
    ctypeRNum t14 = -t10;
    ctypeRNum t15 = t6 + t10;
    ctypeRNum t16 = t7 + t14;
    ctypeRNum t17 = t11 * t12 * t15;
    ctypeRNum t18 = atan(t17);
    ctypeRNum t19 = t11 * t12 * t16;
    ctypeRNum t20 = atan(t19);
    ctypeRNum t21 = -t18;
    ctypeRNum t25 = pSys_[5] * t5 * 3.141592653589793 * (t18 - x[5]) * (-1.0 / 2.0);
    ctypeRNum t22 = t21 + x[5];
    ctypeRNum t23 = (pSys_[6] * t5 * t20 * 3.141592653589793) / 2.0;
    ctypeRNum t26 = atan(t25);
    ctypeRNum t24 = atan(t23);
    out[0] = x[1];
    out[1] = (p[3] - p[0] * pSys_[4] * F_r_z_ * t9 * t13 * t24 * 2.0 + p[0] * pSys_[3] * F_f_z_ * t3 * t9 * t13 * t26 * 2.0) / pSys_[2];
    out[2] = -x[1] - (t11 * (t4 * (p[1] - p[0] * F_f_z_ * t9 * t13 * t26 * sin(x[5]) * 2.0) - t2 * (p[2] + p[0] * F_r_z_ * t9 * t13 * t24 * 2.0 + p[0] * F_f_z_ * t3 * t9 * t13 * t26 * 2.0))) / pSys_[0];
    out[3] = pSys_[7] * cos(t8);
    out[4] = pSys_[7] * sin(t8);
    out[5] = u[0];
}

/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void VehicleProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{

    ctypeRNum t2 = cos(x[2]);
    ctypeRNum t3 = cos(x[5]);
    ctypeRNum t4 = sin(x[2]);
    ctypeRNum t5 = sin(x[5]);
    ctypeRNum t6 = pSys_[3] + pSys_[4];
    ctypeRNum t7 = pSys_[3] * x[1];
    ctypeRNum t8 = pSys_[4] * x[1];
    ctypeRNum t9 = x[0] + x[2];
    ctypeRNum t10 = 3.141592653589793 * 3.141592653589793;
    ctypeRNum t11 = pSys_[3] * pSys_[3];
    ctypeRNum t12 = pSys_[4] * pSys_[4];
    ctypeRNum t13 = pSys_[5] * pSys_[5];
    ctypeRNum t14 = pSys_[6] * pSys_[6];
    ctypeRNum t15 = pSys_[7] * pSys_[7];
    ctypeRNum t16 = x[1] * x[1];
    ctypeRNum t17 = 1.0 / 3.141592653589793;
    ctypeRNum t21 = 1.0 / pSys_[0];
    ctypeRNum t22 = 1.0 / pSys_[2];
    ctypeRNum t23 = 1.0 / pSys_[7];
    ctypeRNum t18 = cos(t9);
    ctypeRNum t19 = pSys_[7] * t4;
    ctypeRNum t20 = sin(t9);
    ctypeRNum t24 = t4 * t7;
    ctypeRNum t25 = t4 * t8;
    ctypeRNum t26 = 1.0 / t2;
    ctypeRNum t27 = t6 * t6;
    ctypeRNum t28 = 1.0 / t6;
    ctypeRNum t31 = t7 * t7;
    ctypeRNum t32 = t8 * t8;
    ctypeRNum t29 = -t19;
    ctypeRNum t30 = pSys_[7] + t24;
    ctypeRNum t33 = pSys_[7] * t20 * vec[3];
    ctypeRNum t34 = -t25;
    ctypeRNum t35 = t7 * t19 * 2.0;
    ctypeRNum t36 = t8 * t19 * 2.0;
    ctypeRNum t37 = pSys_[7] * t18 * vec[4];
    ctypeRNum t38 = t7 + t19;
    ctypeRNum t39 = -t36;
    ctypeRNum t40 = pSys_[7] + t34;
    ctypeRNum t41 = -t33;
    ctypeRNum t42 = t8 + t29;
    ctypeRNum t43 = t23 * t26 * t38;
    ctypeRNum t44 = t15 + t31 + t35;
    ctypeRNum t45 = atan(t43);
    ctypeRNum t46 = t15 + t32 + t39;
    ctypeRNum t47 = t23 * t26 * t42;
    ctypeRNum t49 = 1.0 / t44;
    ctypeRNum t48 = atan(t47);
    ctypeRNum t50 = 1.0 / t46;
    ctypeRNum t51 = -t45;
    ctypeRNum t54 = (t45 - x[5]) * (t45 - x[5]);
    ctypeRNum t55 = pSys_[5] * t6 * 3.141592653589793 * (t45 - x[5]) * (-1.0 / 2.0);
    ctypeRNum t52 = t48 * t48;
    ctypeRNum t53 = t51 + x[5];
    ctypeRNum t56 = atan(t55);
    ctypeRNum t59 = t10 * t13 * t27 * t54;
    ctypeRNum t57 = t10 * t14 * t27 * t52;
    ctypeRNum t60 = t59 + 4.0;
    ctypeRNum t63 = p[0] * F_f_z_ * t3 * t17 * t28 * t56 * 2.0;
    ctypeRNum t64 = p[0] * F_f_z_ * t5 * t17 * t28 * t56 * 2.0;
    ctypeRNum t58 = t57 + 4.0;
    ctypeRNum t62 = 1.0 / t60;
    ctypeRNum t65 = -t64;
    ctypeRNum t61 = 1.0 / t58;
    out[0] = t37 + t41;
    out[1] = vec[0] - vec[1] * (p[0] * pSys_[6] * pSys_[7] * F_r_z_ * t2 * t12 * t22 * t50 * t61 * 4.0 + p[0] * pSys_[5] * pSys_[7] * F_f_z_ * t2 * t3 * t11 * t22 * t49 * t62 * 4.0) + vec[2] * (t21 * t23 * (t2 * (p[0] * pSys_[4] * pSys_[6] * pSys_[7] * F_r_z_ * t2 * t50 * t61 * 4.0 - p[0] * pSys_[3] * pSys_[5] * pSys_[7] * F_f_z_ * t2 * t3 * t49 * t62 * 4.0) - p[0] * pSys_[3] * pSys_[5] * F_f_z_ * t2 * t5 * t19 * t49 * t62 * 4.0) - 1.0);
    out[2] = t37 + t41 + vec[1] * (p[0] * pSys_[4] * pSys_[6] * pSys_[7] * F_r_z_ * t22 * t40 * t50 * t61 * 4.0 - p[0] * pSys_[3] * pSys_[5] * pSys_[7] * F_f_z_ * t3 * t22 * t30 * t49 * t62 * 4.0) - vec[2] * (t2 * t21 * t23 * (p[0] * pSys_[6] * pSys_[7] * F_r_z_ * t40 * t50 * t61 * 4.0 + p[0] * pSys_[5] * pSys_[7] * F_f_z_ * t3 * t30 * t49 * t62 * 4.0) + t4 * t21 * t23 * (p[2] + t63 + p[0] * F_r_z_ * t17 * t28 * atan((pSys_[6] * t6 * t48 * 3.141592653589793) / 2.0) * 2.0) + t2 * t21 * t23 * (p[1] + t65) + p[0] * pSys_[5] * F_f_z_ * t4 * t5 * t21 * t30 * t49 * t62 * 4.0);
    out[3] = 0.0;
    out[4] = 0.0;
    out[5] = -vec[2] * (t2 * t21 * t23 * (t64 - p[0] * pSys_[5] * F_f_z_ * t3 * t62 * 4.0) - t4 * t21 * t23 * (t63 + p[0] * pSys_[5] * F_f_z_ * t5 * t62 * 4.0)) + vec[1] * (p[0] * pSys_[3] * pSys_[5] * F_f_z_ * t3 * t22 * t62 * 4.0 - p[0] * pSys_[3] * F_f_z_ * t5 * t17 * t22 * t28 * t56 * 2.0);
}

/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void VehicleProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = vec[5];
}

/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void VehicleProblemDescription::dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{

    ctypeRNum t2 = cos(x[2]);
    ctypeRNum t3 = sin(x[2]);
    ctypeRNum t4 = pSys_[3] + pSys_[4];
    ctypeRNum t5 = pSys_[3] * x[1];
    ctypeRNum t6 = pSys_[4] * x[1];
    ctypeRNum t7 = 1.0 / 3.141592653589793;
    ctypeRNum t9 = 1.0 / pSys_[0];
    ctypeRNum t10 = 1.0 / pSys_[2];
    ctypeRNum t11 = 1.0 / pSys_[7];
    ctypeRNum t8 = pSys_[7] * t3;
    ctypeRNum t12 = 1.0 / t2;
    ctypeRNum t13 = 1.0 / t4;
    ctypeRNum t14 = -t8;
    ctypeRNum t15 = t5 + t8;
    ctypeRNum t16 = t6 + t14;
    ctypeRNum t17 = t11 * t12 * t15;
    ctypeRNum t18 = atan(t17);
    ctypeRNum t19 = t11 * t12 * t16;
    ctypeRNum t20 = atan(t19);
    ctypeRNum t21 = -t18;
    ctypeRNum t25 = pSys_[5] * t4 * 3.141592653589793 * (t18 - x[5]) * (-1.0 / 2.0);
    ctypeRNum t22 = t21 + x[5];
    ctypeRNum t23 = (pSys_[6] * t4 * t20 * 3.141592653589793) / 2.0;
    ctypeRNum t26 = atan(t25);
    ctypeRNum t24 = atan(t23);
    out[0] = -t7 * t9 * t10 * t13 * (pSys_[0] * pSys_[4] * F_r_z_ * t24 * vec[1] * 2.0 - pSys_[0] * pSys_[3] * F_f_z_ * t26 * vec[1] * cos(x[5]) * 2.0) + t7 * t9 * t10 * t11 * t13 * (pSys_[2] * F_f_z_ * t26 * vec[2] * cos(x[2] - x[5]) * 2.0 + pSys_[2] * F_r_z_ * t2 * t24 * vec[2] * 2.0);
    out[1] = -t3 * t9 * t11 * vec[2];
    out[2] = t2 * t9 * t11 * vec[2];
    out[3] = t10 * vec[1];
}

/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
    -------------------------------------------------- **/
void VehicleProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = pCost_[12] * (u[0] - udes[0]) * (u[0] - udes[0]) + pCost_[0] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[1] * (x[1] - xdes[1]) * (x[1] - xdes[1]) + pCost_[2] * (x[2] - xdes[2]) * (x[2] - xdes[2]) + pCost_[3] * (x[3] - xdes[3]) * (x[3] - xdes[3]) + pCost_[4] * (x[4] - xdes[4]) * (x[4] - xdes[4]) + pCost_[5] * (x[5] - xdes[5]) * (x[5] - xdes[5]);
}

/** Gradient dl/dx **/
void VehicleProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = pCost_[0] * (x[0] - xdes[0]) * 2.0;
    out[1] = pCost_[1] * (x[1] - xdes[1]) * 2.0;
    out[2] = pCost_[2] * (x[2] - xdes[2]) * 2.0;
    out[3] = pCost_[3] * (x[3] - xdes[3]) * 2.0;
    out[4] = pCost_[4] * (x[4] - xdes[4]) * 2.0;
    out[5] = pCost_[5] * (x[5] - xdes[5]) * 2.0;
}

/** Gradient dl/du **/
void VehicleProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = pCost_[12] * (u[0] - udes[0]) * 2.0;
}

/** Gradient dl/dp **/
void VehicleProblemDescription::dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
}

/** Terminal cost V(T,x(T),p,xdes,userparam)
    ---------------------------------------- **/
void VehicleProblemDescription::Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = pCost_[6] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[7] * (x[1] - xdes[1]) * (x[1] - xdes[1]) + pCost_[8] * (x[2] - xdes[2]) * (x[2] - xdes[2]) + pCost_[9] * (x[3] - xdes[3]) * (x[3] - xdes[3]) + pCost_[10] * (x[4] - xdes[4]) * (x[4] - xdes[4]) + pCost_[11] * (x[5] - xdes[5]) * (x[5] - xdes[5]);
}

/** Gradient dV/dx **/
void VehicleProblemDescription::dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = pCost_[6] * (x[0] - xdes[0]) * 2.0;
    out[1] = pCost_[7] * (x[1] - xdes[1]) * 2.0;
    out[2] = pCost_[8] * (x[2] - xdes[2]) * 2.0;
    out[3] = pCost_[9] * (x[3] - xdes[3]) * 2.0;
    out[4] = pCost_[10] * (x[4] - xdes[4]) * 2.0;
    out[5] = pCost_[11] * (x[5] - xdes[5]) * 2.0;
}

/** Gradient dV/dp **/
void VehicleProblemDescription::dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
}

/** Gradient dV/dT **/
void VehicleProblemDescription::dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = 0.0;
}

/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0
    ------------------------------------------------------ **/
void VehicleProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    ctypeRNum t2 = 0.1;
    ctypeRNum t3 = -t2;
    out[0] = x[4] - pCon_[0];
    out[1] = t3 + x[5];
    out[2] = t3 - x[5];
}

/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void VehicleProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = vec[0];
    out[5] = vec[1] - vec[2];
}

/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void VehicleProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0.0;
}

/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void VehicleProblemDescription::dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0.0;
    out[1] = 0.0;
    out[2] = 0.0;
    out[3] = 0.0;
}