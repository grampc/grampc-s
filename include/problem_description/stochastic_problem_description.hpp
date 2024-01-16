#ifndef STOCHASTIC_PROBLEM_DESCRIPTION
#define STOCHASTIC_PROBLEM_DESCRIPTION

#include <vector>
#include <memory>
#include <cmath>

extern "C"
{
#include "grampc.h"
}

namespace grampc
{
    class StochasticProblemDescription
    {
    public:
        virtual ~StochasticProblemDescription() {}

        /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
			inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
		virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) = 0;


		/** System function f(t,x,u,p,userparam)
		------------------------------------ **/
		virtual void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) = 0;

		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		virtual void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) = 0;

		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		virtual void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
		virtual void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) {}


		/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
		virtual void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) {}

		/** Gradient dl/dx **/
		virtual void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) {}

		/** Gradient dl/du **/
		virtual void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) {}

		/** Gradient dl/dp **/
		virtual void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) {}


		/** Terminal cost V(T,x,p) */
		virtual void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) {}

		/** Gradient dV/dx **/
		virtual void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) {}

		/** Gradient dV/dp **/
		virtual void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) {}

		/** Gradient dV/dT **/
		virtual void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) {}


		/** Inequality constraints h(t,x,u,p) < 0 
        ----------------------------------------- **/
		virtual void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		virtual void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}


        /** Equality constraints g(t,x,u,p) = 0 
         * ------------------------------------ */
		virtual void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
		virtual void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) {}

        /** Terminal equality constraints gT(T,x,p) = 0 
         * -------------------------------------------- */
		virtual void gTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) {}

		/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
		virtual void dgTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
		virtual void dgTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
		virtual void dgTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}


		/** Terminal inequality constraints hT(T,x,p) < 0 
         * ---------------------------------------------- */
		virtual void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) {}

		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}

		/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
		virtual void dhTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}
        
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) {}


		/** Additional functions required for semi-implicit systems
		M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
		------------------------------------------------------- **/
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian df/dt **/
		virtual void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian d(dH/dx)/dt  **/
		virtual void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *adj, ctypeRNum *p) {}

		/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mfct(typeRNum *out) {}

		/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mtrans(typeRNum *out) {}


        /** Additional functions required for Taylor-SMPC */
        /*------------------------------------------------*/

        /** Jacobian df/dp in vector form  **/
        virtual void dfdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Hessian d(df/dx)/dx in vector form  */
        virtual void dfdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Jacobian d(df/dx)/dp in vector form  */
        virtual void dfdxdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Hessian d(df/dx)/du in vector form  */
        virtual void dfdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Jacobian d(df/dp)/du in vector form  */
        virtual void dfdpdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Jacobian dh/dx in vector form  **/
        virtual void dhdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Jacobian dh/du in vector form  **/
        virtual void dhdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Hessian d(dh/dx)/dx in vector form  **/
        virtual void dhdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

        /** Jacobian d(dh/dx)/du in vector form **/
        virtual void dhdxdu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) {}

		/** Jacobian dhT/dx in vector form  **/
        virtual void dhTdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) {}

		/** Hessian d(dhT/dx)/dx in vector form  **/
        virtual void dhTdxdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) {}

		/** Hessian d(dhT/dx)/dT in vector form  **/
        virtual void dhTdxdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) {}
    };
    
    // Alias
    typedef std::shared_ptr<StochasticProblemDescription> StochasticProblemDescriptionPtr;
    typedef std::shared_ptr<const StochasticProblemDescription> StochasticProblemDescriptionConstPtr;
}


#endif // STOCHASTIC_PROBLEM_DESCRIPTION