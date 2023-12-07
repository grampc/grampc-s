#include "vehicle_problem_description.hpp"
#include <iostream>

//namespace grampc
//{
    VehicleProblemDescription::VehicleProblemDescription(const std::vector<typeRNum>& pSys,
                                                        const std::vector<typeRNum>& pCost,
                                                        const std::vector<typeRNum>& pCon)
        : pSys_(pSys), pCost_(pCost), pCon_(pCon)
    {
        F_f_z_ = pSys_[0] * pSys_[1] * pSys_[4] / (pSys_[3] + pSys_[4]);
        F_r_z_ = pSys_[0] * pSys_[1] * pSys_[3] / (pSys_[3] + pSys_[4]);
    }

    /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng),
        inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
    void VehicleProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
    {
        *Nx  = 6;
        *Nu  = 1;
        *Np  = 4;
        *Nh  = 3;
        *Ng  = 0;
        *NgT = 0;
        *NhT = 0;
    }


    /** System function f(t,x,u,p,userparam)
        ------------------------------------ **/
    void VehicleProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = cos(x[5]);
        ctypeRNum  t4 = sin(x[2]);
        ctypeRNum  t5 = pSys_[3]+pSys_[4];
        ctypeRNum  t6 = pSys_[3]*x[1];
        ctypeRNum  t7 = pSys_[4]*x[1];
        ctypeRNum  t8 = x[0]+x[2];
        ctypeRNum  t9 = 1.0/3.141592653589793;
        ctypeRNum  t11 = 1.0/pSys_[7];
        ctypeRNum  t10 = pSys_[7]*t4;
        ctypeRNum  t12 = 1.0/t2;
        ctypeRNum  t13 = 1.0/t5;
        ctypeRNum  t14 = -t10;
        ctypeRNum  t15 = t6+t10;
        ctypeRNum  t16 = t7+t14;
        ctypeRNum  t17 = t11*t12*t15;
        ctypeRNum  t18 = atan(t17);
        ctypeRNum  t19 = t11*t12*t16;
        ctypeRNum  t20 = atan(t19);
        ctypeRNum  t21 = -t18;
        ctypeRNum  t25 = pSys_[5]*t5*3.141592653589793*(t18-x[5])*(-1.0/2.0);
        ctypeRNum  t22 = t21+x[5];
        ctypeRNum  t23 = (pSys_[6]*t5*t20*3.141592653589793)/2.0;
        ctypeRNum  t26 = atan(t25);
        ctypeRNum  t24 = atan(t23);
        out[0] = x[1];
        out[1] = (p[3]-p[0]*pSys_[4]*F_r_z_*t9*t13*t24*2.0+p[0]*pSys_[3]*F_f_z_*t3*t9*t13*t26*2.0)/pSys_[2];
        out[2] = -x[1]-(t11*(t4*(p[1]-p[0]*F_f_z_*t9*t13*t26*sin(x[5])*2.0)-t2*(p[2]+p[0]*F_r_z_*t9*t13*t24*2.0+p[0]*F_f_z_*t3*t9*t13*t26*2.0)))/pSys_[0];
        out[3] = pSys_[7]*cos(t8);
        out[4] = pSys_[7]*sin(t8);
        out[5] = u[0];
    }


    /** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
    void VehicleProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = cos(x[5]);
        ctypeRNum  t4 = sin(x[2]);
        ctypeRNum  t5 = sin(x[5]);
        ctypeRNum  t6 = pSys_[3]+pSys_[4];
        ctypeRNum  t7 = pSys_[3]*x[1];
        ctypeRNum  t8 = pSys_[4]*x[1];
        ctypeRNum  t9 = x[0]+x[2];
        ctypeRNum  t10 = 3.141592653589793*3.141592653589793;
        ctypeRNum  t11 = pSys_[3]*pSys_[3];
        ctypeRNum  t12 = pSys_[4]*pSys_[4];
        ctypeRNum  t13 = pSys_[5]*pSys_[5];
        ctypeRNum  t14 = pSys_[6]*pSys_[6];
        ctypeRNum  t15 = pSys_[7]*pSys_[7];
        ctypeRNum  t16 = x[1]*x[1];
        ctypeRNum  t17 = 1.0/3.141592653589793;
        ctypeRNum  t21 = 1.0/pSys_[0];
        ctypeRNum  t22 = 1.0/pSys_[2];
        ctypeRNum  t23 = 1.0/pSys_[7];
        ctypeRNum  t18 = cos(t9);
        ctypeRNum  t19 = pSys_[7]*t4;
        ctypeRNum  t20 = sin(t9);
        ctypeRNum  t24 = t4*t7;
        ctypeRNum  t25 = t4*t8;
        ctypeRNum  t26 = 1.0/t2;
        ctypeRNum  t27 = t6*t6;
        ctypeRNum  t28 = 1.0/t6;
        ctypeRNum  t31 = t7*t7;
        ctypeRNum  t32 = t8*t8;
        ctypeRNum  t29 = -t19;
        ctypeRNum  t30 = pSys_[7]+t24;
        ctypeRNum  t33 = pSys_[7]*t20*vec[3];
        ctypeRNum  t34 = -t25;
        ctypeRNum  t35 = t7*t19*2.0;
        ctypeRNum  t36 = t8*t19*2.0;
        ctypeRNum  t37 = pSys_[7]*t18*vec[4];
        ctypeRNum  t38 = t7+t19;
        ctypeRNum  t39 = -t36;
        ctypeRNum  t40 = pSys_[7]+t34;
        ctypeRNum  t41 = -t33;
        ctypeRNum  t42 = t8+t29;
        ctypeRNum  t43 = t23*t26*t38;
        ctypeRNum  t44 = t15+t31+t35;
        ctypeRNum  t45 = atan(t43);
        ctypeRNum  t46 = t15+t32+t39;
        ctypeRNum  t47 = t23*t26*t42;
        ctypeRNum  t49 = 1.0/t44;
        ctypeRNum  t48 = atan(t47);
        ctypeRNum  t50 = 1.0/t46;
        ctypeRNum  t51 = -t45;
        ctypeRNum  t54 = (t45-x[5])*(t45-x[5]);
        ctypeRNum  t55 = pSys_[5]*t6*3.141592653589793*(t45-x[5])*(-1.0/2.0);
        ctypeRNum  t52 = t48*t48;
        ctypeRNum  t53 = t51+x[5];
        ctypeRNum  t56 = atan(t55);
        ctypeRNum  t59 = t10*t13*t27*t54;
        ctypeRNum  t57 = t10*t14*t27*t52;
        ctypeRNum  t60 = t59+4.0;
        ctypeRNum  t63 = p[0]*F_f_z_*t3*t17*t28*t56*2.0;
        ctypeRNum  t64 = p[0]*F_f_z_*t5*t17*t28*t56*2.0;
        ctypeRNum  t58 = t57+4.0;
        ctypeRNum  t62 = 1.0/t60;
        ctypeRNum  t65 = -t64;
        ctypeRNum  t61 = 1.0/t58;
        out[0] = t37+t41;
        out[1] = vec[0]-vec[1]*(p[0]*pSys_[6]*pSys_[7]*F_r_z_*t2*t12*t22*t50*t61*4.0+p[0]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t11*t22*t49*t62*4.0)+vec[2]*(t21*t23*(t2*(p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t2*t50*t61*4.0-p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t49*t62*4.0)-p[0]*pSys_[3]*pSys_[5]*F_f_z_*t2*t5*t19*t49*t62*4.0)-1.0);
        out[2] = t37+t41+vec[1]*(p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t22*t40*t50*t61*4.0-p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t3*t22*t30*t49*t62*4.0)-vec[2]*(t2*t21*t23*(p[0]*pSys_[6]*pSys_[7]*F_r_z_*t40*t50*t61*4.0+p[0]*pSys_[5]*pSys_[7]*F_f_z_*t3*t30*t49*t62*4.0)+t4*t21*t23*(p[2]+t63+p[0]*F_r_z_*t17*t28*atan((pSys_[6]*t6*t48*3.141592653589793)/2.0)*2.0)+t2*t21*t23*(p[1]+t65)+p[0]*pSys_[5]*F_f_z_*t4*t5*t21*t30*t49*t62*4.0);
        out[3] = 0.0;
        out[4] = 0.0;
        out[5] = -vec[2]*(t2*t21*t23*(t64-p[0]*pSys_[5]*F_f_z_*t3*t62*4.0)-t4*t21*t23*(t63+p[0]*pSys_[5]*F_f_z_*t5*t62*4.0))+vec[1]*(p[0]*pSys_[3]*pSys_[5]*F_f_z_*t3*t22*t62*4.0-p[0]*pSys_[3]*F_f_z_*t5*t17*t22*t28*t56*2.0);
    }


    /** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
    void VehicleProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {
        out[0] = vec[5];
    }


    /** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
    void VehicleProblemDescription::dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = sin(x[2]);
        ctypeRNum  t4 = pSys_[3]+pSys_[4];
        ctypeRNum  t5 = pSys_[3]*x[1];
        ctypeRNum  t6 = pSys_[4]*x[1];
        ctypeRNum  t7 = 1.0/3.141592653589793;
        ctypeRNum  t9 = 1.0/pSys_[0];
        ctypeRNum  t10 = 1.0/pSys_[2];
        ctypeRNum  t11 = 1.0/pSys_[7];
        ctypeRNum  t8 = pSys_[7]*t3;
        ctypeRNum  t12 = 1.0/t2;
        ctypeRNum  t13 = 1.0/t4;
        ctypeRNum  t14 = -t8;
        ctypeRNum  t15 = t5+t8;
        ctypeRNum  t16 = t6+t14;
        ctypeRNum  t17 = t11*t12*t15;
        ctypeRNum  t18 = atan(t17);
        ctypeRNum  t19 = t11*t12*t16;
        ctypeRNum  t20 = atan(t19);
        ctypeRNum  t21 = -t18;
        ctypeRNum  t25 = pSys_[5]*t4*3.141592653589793*(t18-x[5])*(-1.0/2.0);
        ctypeRNum  t22 = t21+x[5];
        ctypeRNum  t23 = (pSys_[6]*t4*t20*3.141592653589793)/2.0;
        ctypeRNum  t26 = atan(t25);
        ctypeRNum  t24 = atan(t23);
        out[0] = -t7*t9*t10*t13*(pSys_[0]*pSys_[4]*F_r_z_*t24*vec[1]*2.0-pSys_[0]*pSys_[3]*F_f_z_*t26*vec[1]*cos(x[5])*2.0)+t7*t9*t10*t11*t13*(pSys_[2]*F_f_z_*t26*vec[2]*cos(x[2]-x[5])*2.0+pSys_[2]*F_r_z_*t2*t24*vec[2]*2.0);
        out[1] = -t3*t9*t11*vec[2];
        out[2] = t2*t9*t11*vec[2];
        out[3] = t10*vec[1];
    }


    /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
        -------------------------------------------------- **/
    void VehicleProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        out[0] = pCost_[12]*(u[0]-udes[0])*(u[0]-udes[0])+pCost_[0]*(x[0]-xdes[0])*(x[0]-xdes[0])+pCost_[1]*(x[1]-xdes[1])*(x[1]-xdes[1])+pCost_[2]*(x[2]-xdes[2])*(x[2]-xdes[2])+pCost_[3]*(x[3]-xdes[3])*(x[3]-xdes[3])+pCost_[4]*(x[4]-xdes[4])*(x[4]-xdes[4])+pCost_[5]*(x[5]-xdes[5])*(x[5]-xdes[5]);
    }


    /** Gradient dl/dx **/
    void VehicleProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        out[0] = pCost_[0]*(x[0]-xdes[0])*2.0;
        out[1] = pCost_[1]*(x[1]-xdes[1])*2.0;
        out[2] = pCost_[2]*(x[2]-xdes[2])*2.0;
        out[3] = pCost_[3]*(x[3]-xdes[3])*2.0;
        out[4] = pCost_[4]*(x[4]-xdes[4])*2.0;
        out[5] = pCost_[5]*(x[5]-xdes[5])*2.0;
    }


    /** Gradient dl/du **/
    void VehicleProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
    {
        out[0] = pCost_[12]*(u[0]-udes[0])*2.0;
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
        out[0] = pCost_[6]*(x[0]-xdes[0])*(x[0]-xdes[0])+pCost_[7]*(x[1]-xdes[1])*(x[1]-xdes[1])+pCost_[8]*(x[2]-xdes[2])*(x[2]-xdes[2])+pCost_[9]*(x[3]-xdes[3])*(x[3]-xdes[3])+pCost_[10]*(x[4]-xdes[4])*(x[4]-xdes[4])+pCost_[11]*(x[5]-xdes[5])*(x[5]-xdes[5]);
    }


    /** Gradient dV/dx **/
    void VehicleProblemDescription::dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
    {
        out[0] = pCost_[6]*(x[0]-xdes[0])*2.0;
        out[1] = pCost_[7]*(x[1]-xdes[1])*2.0;
        out[2] = pCost_[8]*(x[2]-xdes[2])*2.0;
        out[3] = pCost_[9]*(x[3]-xdes[3])*2.0;
        out[4] = pCost_[10]*(x[4]-xdes[4])*2.0;
        out[5] = pCost_[11]*(x[5]-xdes[5])*2.0;
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
        ctypeRNum  t2 = 0.1;
        ctypeRNum  t3 = -t2;
        out[0] = x[4]-pCon_[0];
        out[1] = t3+x[5];
        out[2] = t3-x[5];
    }


    /** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
    void VehicleProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
    {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        out[3] = 0.0;
        out[4] = vec[0];
        out[5] = vec[1]-vec[2];
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

    /** Jacobian df/dx in vector form (column-wise) **/
    void VehicleProblemDescription::dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = cos(x[5]);
        ctypeRNum  t4 = sin(x[2]);
        ctypeRNum  t5 = sin(x[5]);
        ctypeRNum  t6 = pSys_[3]+pSys_[4];
        ctypeRNum  t7 = pSys_[3]*x[1];
        ctypeRNum  t8 = pSys_[4]*x[1];
        ctypeRNum  t9 = x[0]+x[2];
        ctypeRNum  t10 = 3.141592653589793*3.141592653589793;
        ctypeRNum  t11 = pSys_[3]*pSys_[3];
        ctypeRNum  t12 = pSys_[4]*pSys_[4];
        ctypeRNum  t13 = pSys_[5]*pSys_[5];
        ctypeRNum  t14 = pSys_[6]*pSys_[6];
        ctypeRNum  t15 = pSys_[7]*pSys_[7];
        ctypeRNum  t16 = x[1]*x[1];
        ctypeRNum  t17 = 1.0/3.141592653589793;
        ctypeRNum  t21 = 1.0/pSys_[0];
        ctypeRNum  t22 = 1.0/pSys_[2];
        ctypeRNum  t23 = 1.0/pSys_[7];
        ctypeRNum  t18 = cos(t9);
        ctypeRNum  t19 = pSys_[7]*t4;
        ctypeRNum  t20 = sin(t9);
        ctypeRNum  t24 = t4*t7;
        ctypeRNum  t25 = t4*t8;
        ctypeRNum  t26 = 1.0/t2;
        ctypeRNum  t27 = t6*t6;
        ctypeRNum  t30 = 1.0/t6;
        ctypeRNum  t33 = t7*t7;
        ctypeRNum  t34 = t8*t8;
        ctypeRNum  t28 = pSys_[7]*t18;
        ctypeRNum  t29 = pSys_[7]*t20;
        ctypeRNum  t31 = -t19;
        ctypeRNum  t32 = pSys_[7]+t24;
        ctypeRNum  t35 = -t25;
        ctypeRNum  t36 = t7*t19*2.0;
        ctypeRNum  t37 = t8*t19*2.0;
        ctypeRNum  t38 = t7+t19;
        ctypeRNum  t39 = -t29;
        ctypeRNum  t40 = -t37;
        ctypeRNum  t41 = pSys_[7]+t35;
        ctypeRNum  t42 = t8+t31;
        ctypeRNum  t43 = t23*t26*t38;
        ctypeRNum  t44 = t15+t33+t36;
        ctypeRNum  t45 = atan(t43);
        ctypeRNum  t46 = t15+t34+t40;
        ctypeRNum  t47 = t23*t26*t42;
        ctypeRNum  t49 = 1.0/t44;
        ctypeRNum  t48 = atan(t47);
        ctypeRNum  t50 = 1.0/t46;
        ctypeRNum  t51 = -t45;
        ctypeRNum  t54 = (t45-x[5])*(t45-x[5]);
        ctypeRNum  t55 = pSys_[5]*t6*3.141592653589793*(t45-x[5])*(-1.0/2.0);
        ctypeRNum  t52 = t48*t48;
        ctypeRNum  t53 = t51+x[5];
        ctypeRNum  t56 = atan(t55);
        ctypeRNum  t59 = t10*t13*t27*t54;
        ctypeRNum  t57 = t10*t14*t27*t52;
        ctypeRNum  t60 = t59+4.0;
        ctypeRNum  t63 = p[0]*F_f_z_*t3*t17*t30*t56*2.0;
        ctypeRNum  t64 = p[0]*F_f_z_*t5*t17*t30*t56*2.0;
        ctypeRNum  t58 = t57+4.0;
        ctypeRNum  t62 = 1.0/t60;
        ctypeRNum  t65 = -t64;
        ctypeRNum  t61 = 1.0/t58;
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        out[3] = t39;
        out[4] = t28;
        out[5] = 0.0;
        out[6] = 1.0;
        out[7] = p[0]*pSys_[6]*pSys_[7]*F_r_z_*t2*t12*t22*t50*t61*-4.0-p[0]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t11*t22*t49*t62*4.0;
        out[8] = t21*t23*(t2*(p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t2*t50*t61*4.0-p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t49*t62*4.0)-p[0]*pSys_[3]*pSys_[5]*F_f_z_*t2*t5*t19*t49*t62*4.0)-1.0;
        out[9] = 0.0;
        out[10] = 0.0;
        out[11] = 0.0;
        out[12] = 0.0;
        out[13] = p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t22*t41*t50*t61*4.0-p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t3*t22*t32*t49*t62*4.0;
        out[14] = -t2*t21*t23*(p[0]*pSys_[6]*pSys_[7]*F_r_z_*t41*t50*t61*4.0+p[0]*pSys_[5]*pSys_[7]*F_f_z_*t3*t32*t49*t62*4.0)-t4*t21*t23*(p[2]+t63+p[0]*F_r_z_*t17*t30*atan((pSys_[6]*t6*t48*3.141592653589793)/2.0)*2.0)-t2*t21*t23*(p[1]+t65)-p[0]*pSys_[5]*F_f_z_*t4*t5*t21*t32*t49*t62*4.0;
        out[15] = t39;
        out[16] = t28;
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
        out[31] = p[0]*pSys_[3]*pSys_[5]*F_f_z_*t3*t22*t62*4.0-p[0]*pSys_[3]*F_f_z_*t5*t17*t22*t30*t56*2.0;
        out[32] = -t2*t21*t23*(t64-p[0]*pSys_[5]*F_f_z_*t3*t62*4.0)+t4*t21*t23*(t63+p[0]*pSys_[5]*F_f_z_*t5*t62*4.0);
        out[33] = 0.0;
        out[34] = 0.0;
        out[35] = 0.0;
    }


    /** Hessian d(df/dx)/dx in vector form (column-wise) **/
    void VehicleProblemDescription::dfdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = cos(x[5]);
        ctypeRNum  t4 = sin(x[2]);
        ctypeRNum  t5 = sin(x[5]);
        ctypeRNum  t6 = pSys_[3]+pSys_[4];
        ctypeRNum  t7 = pSys_[3]*x[1];
        ctypeRNum  t8 = pSys_[4]*x[1];
        ctypeRNum  t9 = x[0]+x[2];
        ctypeRNum  t10 = 3.141592653589793*3.141592653589793;
        ctypeRNum  t11 = pSys_[3]*pSys_[3];
        ctypeRNum  t12 = pSys_[3]*pSys_[3]*pSys_[3];
        ctypeRNum  t13 = pSys_[4]*pSys_[4];
        ctypeRNum  t14 = pSys_[4]*pSys_[4]*pSys_[4];
        ctypeRNum  t15 = pSys_[5]*pSys_[5];
        ctypeRNum  t16 = pSys_[5]*pSys_[5]*pSys_[5];
        ctypeRNum  t17 = pSys_[6]*pSys_[6];
        ctypeRNum  t18 = pSys_[6]*pSys_[6]*pSys_[6];
        ctypeRNum  t19 = pSys_[7]*pSys_[7];
        ctypeRNum  t20 = x[5]*2.0;
        ctypeRNum  t21 = x[1]*x[1];
        ctypeRNum  t23 = 1.0/3.141592653589793;
        ctypeRNum  t27 = 1.0/pSys_[0];
        ctypeRNum  t28 = 1.0/pSys_[2];
        ctypeRNum  t29 = 1.0/pSys_[7];
        ctypeRNum  t22 = t2*t2;
        ctypeRNum  t24 = cos(t9);
        ctypeRNum  t25 = pSys_[7]*t4;
        ctypeRNum  t26 = sin(t9);
        ctypeRNum  t30 = t4*t7;
        ctypeRNum  t31 = t4*t8;
        ctypeRNum  t32 = 1.0/t2;
        ctypeRNum  t33 = t6*t6;
        ctypeRNum  t36 = 1.0/t6;
        ctypeRNum  t39 = t7*t7;
        ctypeRNum  t40 = t8*t8;
        ctypeRNum  t34 = pSys_[7]*t24;
        ctypeRNum  t35 = pSys_[7]*t26;
        ctypeRNum  t37 = -t25;
        ctypeRNum  t38 = pSys_[7]+t30;
        ctypeRNum  t41 = -t31;
        ctypeRNum  t42 = t7*t25*2.0;
        ctypeRNum  t43 = t8*t25*2.0;
        ctypeRNum  t44 = t7+t25;
        ctypeRNum  t45 = -t34;
        ctypeRNum  t46 = -t35;
        ctypeRNum  t47 = -t43;
        ctypeRNum  t48 = pSys_[7]+t41;
        ctypeRNum  t49 = t38*t38;
        ctypeRNum  t50 = t8+t37;
        ctypeRNum  t52 = t29*t32*t44;
        ctypeRNum  t53 = t19+t39+t42;
        ctypeRNum  t51 = t48*t48;
        ctypeRNum  t54 = atan(t52);
        ctypeRNum  t55 = t19+t40+t47;
        ctypeRNum  t56 = t29*t32*t50;
        ctypeRNum  t58 = 1.0/t53;
        ctypeRNum  t57 = atan(t56);
        ctypeRNum  t59 = t58*t58;
        ctypeRNum  t60 = t54*2.0;
        ctypeRNum  t61 = 1.0/t55;
        ctypeRNum  t63 = -t54;
        ctypeRNum  t68 = (t54-x[5])*(t54-x[5]);
        ctypeRNum  t69 = pSys_[5]*t6*3.141592653589793*(t54-x[5])*(-1.0/2.0);
        ctypeRNum  t62 = t61*t61;
        ctypeRNum  t64 = -t60;
        ctypeRNum  t65 = t57*t57;
        ctypeRNum  t66 = t63+x[5];
        ctypeRNum  t70 = atan(t69);
        ctypeRNum  t73 = t10*t15*t33*t68;
        ctypeRNum  t67 = t20+t64;
        ctypeRNum  t71 = t10*t17*t33*t65;
        ctypeRNum  t74 = t73+4.0;
        ctypeRNum  t80 = p[0]*F_f_z_*t3*t23*t36*t70*2.0;
        ctypeRNum  t82 = p[0]*F_f_z_*t5*t23*t36*t70*2.0;
        ctypeRNum  t72 = t71+4.0;
        ctypeRNum  t77 = 1.0/t74;
        ctypeRNum  t83 = -t82;
        ctypeRNum  t75 = 1.0/t72;
        ctypeRNum  t78 = t77*t77;
        ctypeRNum  t79 = p[0]*pSys_[5]*F_f_z_*t3*t77*4.0;
        ctypeRNum  t81 = p[0]*pSys_[5]*F_f_z_*t5*t77*4.0;
        ctypeRNum  t89 = p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t58*t77*-4.0;
        ctypeRNum  t94 = p[0]*pSys_[3]*pSys_[5]*F_f_z_*t5*t22*t27*t58*t77*-4.0;
        ctypeRNum  t76 = t75*t75;
        ctypeRNum  t84 = p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t2*t61*t75*4.0;
        ctypeRNum  t85 = p[0]*pSys_[4]*pSys_[6]*F_r_z_*t25*t61*t75*4.0;
        ctypeRNum  t86 = pSys_[3]*pSys_[7]*t2*t58*t79;
        ctypeRNum  t87 = pSys_[3]*pSys_[7]*t2*t58*t81;
        ctypeRNum  t88 = pSys_[3]*t25*t58*t79;
        ctypeRNum  t90 = p[0]*pSys_[6]*F_r_z_*t13*t25*t28*t61*t75*4.0;
        ctypeRNum  t91 = pSys_[7]*t38*t58*t79;
        ctypeRNum  t92 = pSys_[7]*t38*t58*t81;
        ctypeRNum  t93 = pSys_[3]*t22*t27*t58*t81;
        ctypeRNum  t95 = pSys_[7]*t2*t11*t28*t58*t81;
        ctypeRNum  t96 = t11*t25*t28*t58*t79;
        ctypeRNum  t98 = t80+t81;
        ctypeRNum  t99 = t79+t83;
        ctypeRNum  t104 = p[0]*pSys_[3]*F_f_z_*t2*t3*t10*t16*t19*t33*t38*t59*t78*(t54-x[5])*-8.0;
        ctypeRNum  t106 = p[0]*pSys_[3]*F_f_z_*t2*t3*t10*t16*t19*t33*t38*t59*t78*(t54-x[5])*8.0;
        ctypeRNum  t108 = p[0]*pSys_[3]*F_f_z_*t2*t5*t10*t16*t25*t27*t33*t38*t59*t78*(t54-x[5])*-8.0;
        ctypeRNum  t109 = p[0]*pSys_[3]*F_f_z_*t2*t5*t10*t16*t25*t27*t33*t38*t59*t78*(t54-x[5])*8.0;
        ctypeRNum  t110 = p[0]*F_f_z_*t2*t3*t10*t11*t16*t19*t28*t33*t38*t59*t78*(t54-x[5])*-8.0;
        ctypeRNum  t111 = p[0]*F_f_z_*t2*t3*t10*t11*t16*t19*t28*t33*t38*t59*t78*(t54-x[5])*8.0;
        ctypeRNum  t97 = pSys_[3]*t28*t92;
        ctypeRNum  t100 = t2*t27*t29*t98;
        ctypeRNum  t101 = t4*t27*t29*t99;
        ctypeRNum  t102 = p[0]*pSys_[4]*F_r_z_*t2*t10*t18*t19*t33*t48*t57*t62*t76*8.0;
        ctypeRNum  t105 = p[0]*F_r_z_*t2*t10*t13*t18*t19*t28*t33*t48*t57*t62*t76*8.0;
        ctypeRNum  t112 = t84+t89;
        ctypeRNum  t103 = -t101;
        ctypeRNum  t107 = -t105;
        ctypeRNum  t113 = t4*t27*t29*t112;
        ctypeRNum  t114 = -t113;
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        out[3] = t45;
        out[4] = t46;
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
        out[15] = t45;
        out[16] = t46;
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
        out[43] = p[0]*pSys_[6]*pSys_[7]*F_r_z_*t2*t14*t28*t50*t62*t75*8.0+p[0]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t12*t28*t44*t59*t77*8.0+p[0]*F_r_z_*t10*t14*t18*t19*t22*t28*t33*t57*t62*t76*8.0+p[0]*F_f_z_*t3*t10*t12*t16*t19*t22*t28*t33*t59*t78*(t54-x[5])*8.0;
        out[44] = -t2*t27*t29*(p[0]*pSys_[6]*pSys_[7]*F_r_z_*t2*t13*t50*t62*t75*8.0-p[0]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t11*t44*t59*t77*8.0+p[0]*F_r_z_*t10*t13*t18*t19*t22*t33*t57*t62*t76*8.0-p[0]*F_f_z_*t3*t10*t11*t16*t19*t22*t33*t59*t78*(t54-x[5])*8.0)+p[0]*pSys_[5]*F_f_z_*t2*t4*t5*t11*t27*t44*t59*t77*8.0+p[0]*F_f_z_*t5*t10*t11*t16*t22*t25*t27*t33*t59*t78*(t54-x[5])*8.0;
        out[45] = 0.0;
        out[46] = 0.0;
        out[47] = 0.0;
        out[48] = 0.0;
        out[49] = -t90+t107+t111-p[0]*pSys_[5]*F_f_z_*t3*t11*t25*t28*t58*t77*4.0-p[0]*pSys_[6]*pSys_[7]*F_r_z_*t13*t28*t48*t50*t62*t75*8.0+p[0]*pSys_[5]*pSys_[7]*F_f_z_*t3*t11*t28*t38*t44*t59*t77*8.0;
        out[50] = t94+t109+t114+t2*t27*t29*(t85+t102+t106-p[0]*pSys_[3]*pSys_[5]*F_f_z_*t3*t25*t58*t77*4.0+p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t48*t50*t62*t75*8.0+p[0]*pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t3*t38*t44*t59*t77*8.0)+pSys_[3]*t27*t58*t81*(t22-1.0)+p[0]*pSys_[3]*pSys_[5]*F_f_z_*t4*t5*t27*t38*t44*t59*t77*8.0;
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
        out[67] = t95-p[0]*pSys_[7]*F_f_z_*t2*t3*t10*t11*t16*t28*t33*t58*t78*(t54-x[5])*8.0;
        out[68] = t2*t27*t29*(t87-p[0]*pSys_[3]*pSys_[7]*F_f_z_*t2*t3*t10*t16*t33*t58*t78*(t54-x[5])*8.0)-t4*t27*t29*(t86+p[0]*pSys_[3]*pSys_[7]*F_f_z_*t2*t5*t10*t16*t33*t58*t78*(t54-x[5])*8.0);
        out[69] = 0.0;
        out[70] = 0.0;
        out[71] = 0.0;
        out[72] = 0.0;
        out[73] = 0.0;
        out[74] = 0.0;
        out[75] = t45;
        out[76] = t46;
        out[77] = 0.0;
        out[78] = 0.0;
        out[79] = t90+t96+t107+t111-p[0]*pSys_[6]*F_r_z_*t8*t13*t19*t22*t28*t62*t75*8.0+p[0]*pSys_[5]*F_f_z_*t3*t7*t11*t19*t22*t28*t59*t77*8.0;
        out[80] = t94+t109+t114+t2*t27*t29*(-t85+t88+t102+t106+p[0]*pSys_[4]*pSys_[6]*F_r_z_*t8*t19*t22*t62*t75*8.0+p[0]*pSys_[3]*pSys_[5]*F_f_z_*t3*t7*t19*t22*t59*t77*8.0)+pSys_[3]*(t4*t4)*t27*t58*t81+p[0]*pSys_[3]*pSys_[5]*F_f_z_*t5*t7*t22*t25*t27*t59*t77*8.0;
        out[81] = 0.0;
        out[82] = 0.0;
        out[83] = 0.0;
        out[84] = 0.0;
        out[85] = t7*t28*t89-p[0]*pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t2*t8*t28*t61*t75*4.0+p[0]*pSys_[4]*pSys_[6]*F_r_z_*t2*t8*t19*t28*t48*t62*t75*8.0+p[0]*pSys_[3]*pSys_[5]*F_f_z_*t2*t3*t7*t19*t28*t38*t59*t77*8.0+p[0]*pSys_[4]*F_r_z_*t10*t18*t19*t28*t33*t51*t57*t62*t76*8.0+p[0]*pSys_[3]*F_f_z_*t3*t10*t16*t19*t28*t33*t49*t59*t78*(t54-x[5])*8.0;
        out[86] = t4*t27*t29*(t91+p[0]*pSys_[6]*pSys_[7]*F_r_z_*t48*t61*t75*4.0)*2.0+t2*t27*t29*(p[0]*pSys_[6]*pSys_[7]*F_r_z_*t2*t8*t61*t75*4.0-p[0]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t7*t58*t77*4.0-p[0]*pSys_[6]*F_r_z_*t2*t8*t19*t48*t62*t75*8.0+p[0]*pSys_[5]*F_f_z_*t2*t3*t7*t19*t38*t59*t77*8.0-p[0]*F_r_z_*t10*t18*t19*t33*t51*t57*t62*t76*8.0+p[0]*F_f_z_*t3*t10*t16*t19*t33*t49*t59*t78*(t54-x[5])*8.0)-t2*t27*t29*(p[2]+t80+p[0]*F_r_z_*t23*t36*atan((pSys_[6]*t6*t57*3.141592653589793)/2.0)*2.0)+t4*t27*t29*(p[1]+t83)-p[0]*pSys_[5]*F_f_z_*t2*t5*t27*t30*t58*t77*4.0-p[0]*pSys_[5]*F_f_z_*t2*t5*t27*t38*t58*t77*8.0+p[0]*pSys_[5]*F_f_z_*t2*t5*t7*t25*t27*t38*t59*t77*8.0+p[0]*F_f_z_*t5*t10*t16*t25*t27*t33*t49*t59*t78*(t54-x[5])*8.0;
        out[87] = t45;
        out[88] = t46;
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
        out[103] = t97-p[0]*pSys_[3]*pSys_[7]*F_f_z_*t3*t10*t16*t28*t33*t38*t58*t78*(t54-x[5])*8.0;
        out[104] = t100+t103+t2*t27*t29*(t92-p[0]*pSys_[7]*F_f_z_*t3*t10*t16*t33*t38*t58*t78*(t54-x[5])*8.0)-t4*t27*t29*(t91+p[0]*pSys_[7]*F_f_z_*t5*t10*t16*t33*t38*t58*t78*(t54-x[5])*8.0);
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
        out[187] = t95+p[0]*pSys_[7]*F_f_z_*t2*t3*t10*t11*t16*t28*t33*t58*t67*t78*4.0;
        out[188] = t2*t27*t29*(t87+p[0]*pSys_[3]*pSys_[7]*F_f_z_*t2*t3*t10*t16*t33*t58*t67*t78*4.0)-p[0]*pSys_[3]*pSys_[5]*F_f_z_*t2*t3*t4*t27*t58*t77*4.0+p[0]*pSys_[3]*F_f_z_*t2*t4*t5*t10*t16*t27*t33*t58*t67*t78*4.0;
        out[189] = 0.0;
        out[190] = 0.0;
        out[191] = 0.0;
        out[192] = 0.0;
        out[193] = t97+p[0]*pSys_[3]*pSys_[7]*F_f_z_*t3*t10*t16*t28*t33*t38*t58*t67*t78*4.0;
        out[194] = t100+t103+t2*t27*t29*(t92+p[0]*pSys_[7]*F_f_z_*t3*t10*t16*t33*t38*t58*t67*t78*4.0)-p[0]*pSys_[5]*F_f_z_*t3*t4*t27*t38*t58*t77*4.0+p[0]*F_f_z_*t4*t5*t10*t16*t27*t33*t38*t58*t67*t78*4.0;
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
        out[211] = p[0]*pSys_[3]*pSys_[5]*F_f_z_*t5*t28*t77*-8.0-p[0]*pSys_[3]*F_f_z_*t3*t23*t28*t36*t70*2.0-p[0]*pSys_[3]*F_f_z_*t3*t10*t16*t28*t33*t67*t78*4.0;
        out[212] = -t2*t27*t29*(t80+p[0]*pSys_[5]*F_f_z_*t5*t77*8.0+p[0]*F_f_z_*t3*t10*t16*t33*t67*t78*4.0)-t4*t27*t29*(t82-p[0]*pSys_[5]*F_f_z_*t3*t77*8.0+p[0]*F_f_z_*t5*t10*t16*t33*t67*t78*4.0);
        out[213] = 0.0;
        out[214] = 0.0;
        out[215] = 0.0;
    }


    /** Jacobian d(df/dx)/du in vector form (column-wise **/
    void VehicleProblemDescription::dfdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
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
    }


    /** Jacobian d(df/dx)/dp in vector form (column-wise **/
    void VehicleProblemDescription::dfdxdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = cos(x[5]);
        ctypeRNum  t4 = sin(x[2]);
        ctypeRNum  t5 = sin(x[5]);
        ctypeRNum  t6 = pSys_[3]+pSys_[4];
        ctypeRNum  t7 = pSys_[3]*x[1];
        ctypeRNum  t8 = pSys_[4]*x[1];
        ctypeRNum  t9 = 3.141592653589793*3.141592653589793;
        ctypeRNum  t10 = pSys_[3]*pSys_[3];
        ctypeRNum  t11 = pSys_[4]*pSys_[4];
        ctypeRNum  t12 = pSys_[5]*pSys_[5];
        ctypeRNum  t13 = pSys_[6]*pSys_[6];
        ctypeRNum  t14 = pSys_[7]*pSys_[7];
        ctypeRNum  t15 = x[1]*x[1];
        ctypeRNum  t16 = 1.0/3.141592653589793;
        ctypeRNum  t18 = 1.0/pSys_[0];
        ctypeRNum  t19 = 1.0/pSys_[2];
        ctypeRNum  t20 = 1.0/pSys_[7];
        ctypeRNum  t17 = pSys_[7]*t4;
        ctypeRNum  t21 = t4*t7;
        ctypeRNum  t22 = t4*t8;
        ctypeRNum  t23 = 1.0/t2;
        ctypeRNum  t24 = t6*t6;
        ctypeRNum  t25 = 1.0/t6;
        ctypeRNum  t28 = t7*t7;
        ctypeRNum  t29 = t8*t8;
        ctypeRNum  t26 = -t17;
        ctypeRNum  t27 = pSys_[7]+t21;
        ctypeRNum  t30 = -t22;
        ctypeRNum  t31 = t7*t17*2.0;
        ctypeRNum  t32 = t8*t17*2.0;
        ctypeRNum  t33 = t7+t17;
        ctypeRNum  t34 = -t32;
        ctypeRNum  t35 = pSys_[7]+t30;
        ctypeRNum  t36 = t8+t26;
        ctypeRNum  t37 = t20*t23*t33;
        ctypeRNum  t38 = t14+t28+t31;
        ctypeRNum  t39 = atan(t37);
        ctypeRNum  t40 = t14+t29+t34;
        ctypeRNum  t41 = t20*t23*t36;
        ctypeRNum  t43 = 1.0/t38;
        ctypeRNum  t42 = atan(t41);
        ctypeRNum  t44 = 1.0/t40;
        ctypeRNum  t45 = -t39;
        ctypeRNum  t48 = (t39-x[5])*(t39-x[5]);
        ctypeRNum  t49 = pSys_[5]*t6*3.141592653589793*(t39-x[5])*(-1.0/2.0);
        ctypeRNum  t46 = t42*t42;
        ctypeRNum  t47 = t45+x[5];
        ctypeRNum  t50 = atan(t49);
        ctypeRNum  t53 = t9*t12*t24*t48;
        ctypeRNum  t51 = t9*t13*t24*t46;
        ctypeRNum  t54 = t53+4.0;
        ctypeRNum  t52 = t51+4.0;
        ctypeRNum  t56 = 1.0/t54;
        ctypeRNum  t55 = 1.0/t52;
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
        out[3] = 0.0;
        out[4] = 0.0;
        out[5] = 0.0;
        out[6] = 0.0;
        out[7] = pSys_[6]*pSys_[7]*F_r_z_*t2*t11*t19*t44*t55*-4.0-pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t10*t19*t43*t56*4.0;
        out[8] = t2*t18*t20*(pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t2*t44*t55*4.0-pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t2*t3*t43*t56*4.0)-pSys_[3]*pSys_[5]*F_f_z_*t2*t4*t5*t18*t43*t56*4.0;
        out[9] = 0.0;
        out[10] = 0.0;
        out[11] = 0.0;
        out[12] = 0.0;
        out[13] = pSys_[4]*pSys_[6]*pSys_[7]*F_r_z_*t19*t35*t44*t55*4.0-pSys_[3]*pSys_[5]*pSys_[7]*F_f_z_*t3*t19*t27*t43*t56*4.0;
        out[14] = -t2*t18*t20*(pSys_[6]*pSys_[7]*F_r_z_*t35*t44*t55*4.0+pSys_[5]*pSys_[7]*F_f_z_*t3*t27*t43*t56*4.0)-t4*t16*t18*t20*t25*(F_r_z_*atan((pSys_[6]*t6*t42*3.141592653589793)/2.0)+F_f_z_*t3*t50)*2.0-pSys_[5]*F_f_z_*t4*t5*t18*t27*t43*t56*4.0+F_f_z_*t2*t5*t16*t18*t20*t25*t50*2.0;
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
        out[31] = pSys_[3]*pSys_[5]*F_f_z_*t3*t19*t56*4.0-pSys_[3]*F_f_z_*t5*t16*t19*t25*t50*2.0;
        out[32] = t2*t18*t20*(pSys_[5]*F_f_z_*t3*t56*4.0-F_f_z_*t5*t16*t25*t50*2.0)+t4*t18*t20*(pSys_[5]*F_f_z_*t5*t56*4.0+F_f_z_*t3*t16*t25*t50*2.0);
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
        out[50] = -t2*t18*t20;
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
        out[86] = -t4*t18*t20;
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
    }


    /** Jacobian df/dp in vector form (column-wise **/
    void VehicleProblemDescription::dfdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {


        ctypeRNum  t2 = cos(x[2]);
        ctypeRNum  t3 = sin(x[2]);
        ctypeRNum  t4 = pSys_[3]+pSys_[4];
        ctypeRNum  t5 = pSys_[3]*x[1];
        ctypeRNum  t6 = pSys_[4]*x[1];
        ctypeRNum  t7 = 1.0/3.141592653589793;
        ctypeRNum  t9 = 1.0/pSys_[0];
        ctypeRNum  t10 = 1.0/pSys_[2];
        ctypeRNum  t11 = 1.0/pSys_[7];
        ctypeRNum  t8 = pSys_[7]*t3;
        ctypeRNum  t12 = 1.0/t2;
        ctypeRNum  t13 = 1.0/t4;
        ctypeRNum  t14 = -t8;
        ctypeRNum  t15 = t5+t8;
        ctypeRNum  t16 = t6+t14;
        ctypeRNum  t17 = t11*t12*t15;
        ctypeRNum  t18 = atan(t17);
        ctypeRNum  t19 = t11*t12*t16;
        ctypeRNum  t20 = atan(t19);
        ctypeRNum  t21 = -t18;
        ctypeRNum  t25 = pSys_[5]*t4*3.141592653589793*(t18-x[5])*(-1.0/2.0);
        ctypeRNum  t22 = t21+x[5];
        ctypeRNum  t23 = (pSys_[6]*t4*t20*3.141592653589793)/2.0;
        ctypeRNum  t26 = atan(t25);
        ctypeRNum  t24 = atan(t23);
        out[0] = 0.0;
        out[1] = pSys_[4]*F_r_z_*t7*t10*t13*t24*-2.0+pSys_[3]*F_f_z_*t7*t10*t13*t26*cos(x[5])*2.0;
        out[2] = F_f_z_*t7*t9*t11*t13*t26*cos(x[2]-x[5])*2.0+F_r_z_*t2*t7*t9*t11*t13*t24*2.0;
        out[3] = 0.0;
        out[4] = 0.0;
        out[5] = 0.0;
        out[6] = 0.0;
        out[7] = 0.0;
        out[8] = -t3*t9*t11;
        out[9] = 0.0;
        out[10] = 0.0;
        out[11] = 0.0;
        out[12] = 0.0;
        out[13] = 0.0;
        out[14] = t2*t9*t11;
        out[15] = 0.0;
        out[16] = 0.0;
        out[17] = 0.0;
        out[18] = 0.0;
        out[19] = t10;
        out[20] = 0.0;
        out[21] = 0.0;
        out[22] = 0.0;
        out[23] = 0.0;
    }


    /** Jacobian d(df/dp)/du in vector form (column-wise **/
    void VehicleProblemDescription::dfdpdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
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
    }

    /** Jacobian dh/dx in vector form (column-wise) **/
    void VehicleProblemDescription::dhdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
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
        out[12] = 1.0;
        out[13] = 0.0;
        out[14] = 0.0;
        out[15] = 0.0;
        out[16] = 1.0;
        out[17] = -1.0;
    }


    /** Jacobian dh/du in vector form (column-wise) **/
    void VehicleProblemDescription::dhdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
    {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
    }

    /** Hessian d(dh/dx)/dx in vector form (column-wise) **/
    void VehicleProblemDescription::dhdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
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
    }


    /** Jacobian d(dh/dx)/du in vector form (column-wise) **/
    void VehicleProblemDescription::dhdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
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
    }
//}