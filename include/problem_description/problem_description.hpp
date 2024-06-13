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


#ifndef PROBLEM_DESCRIPTION
#define PROBLEM_DESCRIPTION

#include <vector>
#include <memory>
#include <cmath>

#include "util/grampc_s_constants.hpp"

namespace grampc
{
    class ProblemDescription
    {
    public:
        virtual ~ProblemDescription() {}

        /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
			inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
		virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) = 0;


		/** System function f(t,x,u,p,userparam)
		------------------------------------ **/
		virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) = 0;

		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) = 0;

		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
		virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) {}


		/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
		virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}

		/** Gradient dl/dx **/
		virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}

		/** Gradient dl/du **/
		virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}

		/** Gradient dl/dp **/
		virtual void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}


		/** Terminal cost V(T,x,p) */
		virtual void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}

		/** Gradient dV/dx **/
		virtual void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}

		/** Gradient dV/dp **/
		virtual void dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}

		/** Gradient dV/dT **/
		virtual void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}


		/** Inequality constraints h(t,x,u,p) < 0 
        ----------------------------------------- **/
		virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}


        /** Equality constraints g(t,x,u,p) = 0 
         * ------------------------------------ */
		virtual void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
		virtual void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}

        /** Terminal equality constraints gT(T,x,p) = 0 
         * -------------------------------------------- */
		virtual void gTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}

		/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
		virtual void dgTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
		virtual void dgTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
		virtual void dgTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}


		/** Terminal inequality constraints hT(T,x,p) < 0 
         * ---------------------------------------------- */
		virtual void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}

		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}

		/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
		virtual void dhTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}
        
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}


		/** Additional functions required for semi-implicit systems
		M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
		------------------------------------------------------- **/
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian df/dt **/
		virtual void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian d(dH/dx)/dt  **/
		virtual void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef adj, VectorConstRef p) {}

		/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mfct(VectorRef out) {}

		/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mtrans(VectorRef out) {}


        /** Additional functions required for Taylor-SMPC */
        /*------------------------------------------------*/

        /** Jacobian df/dp in vector form  **/
        virtual void dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Hessian d(df/dx)/dx in vector form  */
        virtual void dfdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Jacobian d(df/dx)/dp in vector form  */
        virtual void dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Hessian d(df/dx)/du in vector form  */
        virtual void dfdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Jacobian d(df/dp)/du in vector form  */
        virtual void dfdpdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Jacobian dh/dx in vector form  **/
        virtual void dhdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Jacobian dh/du in vector form  **/
        virtual void dhdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Hessian d(dh/dx)/dx in vector form  **/
        virtual void dhdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

        /** Jacobian d(dh/dx)/du in vector form **/
        virtual void dhdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}

		/** Jacobian dhT/dx in vector form  **/
        virtual void dhTdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}

		/** Hessian d(dhT/dx)/dx in vector form  **/
        virtual void dhTdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}

		/** Hessian d(dhT/dx)/dT in vector form  **/
        virtual void dhTdxdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}

	public:
		typeInt Nx_;
		typeInt Np_;
		typeInt Nu_;
		typeInt Nh_;
		typeInt NhT_;
    };
    
    // Alias
    typedef std::shared_ptr<ProblemDescription> ProblemDescriptionPtr;
    typedef std::shared_ptr<const ProblemDescription> ProblemDescriptionConstPtr;
}


#endif // PROBLEM_DESCRIPTION