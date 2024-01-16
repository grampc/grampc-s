#ifndef TAYLOR_PROBLEM_DESCRIPTION_HPP
#define TAYLOR_PROBLEM_DESCRIPTION_HPP


#include "problem_description/grampc_interface.hpp"
#include "point_transformation/point_transformation.hpp"
#include "constraint_approx/chance_constraint_approximation.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    // Interface to GRAMPC using a first order Taylor series approximation for uncertainty propagation
    class TaylorProblemDescription : public StochasticProblemDescription
    {
    public:
        // Constructor specifying the type of constraint approximation
        TaylorProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation);

        // Constructor with zero-mean Wiener process and specified covariance matrix
        TaylorProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef covWienerProcess);
        
        /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
			inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;


        /** System function f(t,x,u,p,userparam)
		------------------------------------ **/
        virtual void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        
        // Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx)
        virtual void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) override;
        
        // Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du)
        virtual void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p) override;


        /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
        virtual void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;
        
        // Gradient dl/dx
        virtual void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;
        
        // Gradient dl/du
        virtual void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes) override;


        /** Terminal cost V(T,x,p)
		-------------------------------------------------- **/
        virtual void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;
        
        // Gradient dV/dx
        virtual void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;

        /** Gradient dV/dT **/
		virtual void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes) override;


        /** Inequality constraints h(t,x,u,p) < 0
		-------------------------------------------------- **/
        virtual void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        
        //  Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx)
        virtual void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;

        // Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du)
        virtual void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;


        /** Terminal inequality constraints hT(T,x,p) < 0 
         * ------------------------------------------------ **/
        virtual void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p) override;

		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) override;
        
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec) override;


        // Compute the point-based representation of the initial states and the parameters
        void compute_x0_and_p0(DistributionPtr state, DistributionPtr param);

        // Compute the point-based representation of the initial states and do not replace the parameters
        void compute_x0_and_p0(DistributionPtr state);

        // Return initial state
        ctypeRNum* x0();

        // Return parameters
        ctypeRNum* p0();

    private:  
        typeInt numStates_;
        typeInt numParams_;
        typeInt numControlInputs_;
        typeInt numConstraints_;
        typeInt numTerminalConstraints_;
        
        StochasticProblemDescriptionPtr problemDescription_;
        Vector constraintTighteningCoeff_;
        Matrix covWienerProcess_;
        Matrix x0_;
        Matrix p0_;
        Matrix dfdx_;
        Matrix dfdp_;
        Matrix dfdxdx_;
        Matrix dfdxdp_;
        Matrix dfdxdu_;
        Matrix dfdpdu_;
        Matrix dhdx_;
        Vector dhdvar_;
        Matrix dhdxdx_;
        Matrix dhdxdu_;
        Vector constraintStdDev_;
        Matrix dhTdx_;
        Vector dhTdvar_;
        Matrix dhTdxdx_;
        Matrix dhTdxdT_;
        Vector terminalConstraintStdDev_;

        // Temporary elements
        typeRNum temp_scalar_;
        Vector temp_vec_numStates_;
        Matrix temp_mat_numStates_;
        Matrix temp_mat_numStates_numParams_;
        Matrix temp_mat_numStates_numControlInputs_;
    };
}

#endif // TAYLOR_PROBLEM_DESCRIPTION_HPP
