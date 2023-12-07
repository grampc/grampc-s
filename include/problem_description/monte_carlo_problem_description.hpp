#ifndef MONTE_CARLO_PROBLEM_DESCRIPTION_HPP
#define MONTE_CARLO_PROBLEM_DESCRIPTION_HPP

#include "problem_description/grampc_interface.hpp"
#include "point_transformation/point_transformation.hpp"
#include "constraint_approx/chance_constraint_approximation.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    // Interface to GRAMPC using a Monte-Carlo simulation to represent the uncertainties, which ensures constraint satisfaction for all sampling points
    class MonteCarloProblemDescription : public StochasticProblemDescription
    {
    public:
        // Constructor for a problem description based on a Monte-Carlo simulation 
        MonteCarloProblemDescription(StochasticProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation);
        
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


        /** Inequality constraints h(t,x,u,p) < 0
		-------------------------------------------------- **/
        virtual void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p) override;
        
        //  Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx)
        virtual void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;

        // Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du)
        virtual void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec) override;


        // Compute the point-based representation of the initial states and the parameters
        void compute_x0_and_p0(DistributionPtr state, DistributionPtr param);

        // Return initial state
        ctypeRNum* x0();

        // Return parameters
        ctypeRNum* p0();

    private:  
        typeInt numSigmaPoints_;
        typeInt numStates_;
        typeInt numParams_;
        typeInt numControlInputs_;
        typeInt numConstraints_;
        
        StochasticProblemDescriptionPtr problemDescription_;
        PointTransformationPtr pointTransformation_;
        MultiDistributionPtr stateAndParam_;
        Matrix x0_;
        Matrix p0_;
        RowVector dmean_dPoints_;

        // Temporary vectors
        Vector temp_vec_numInputs_;
        RowVector temp_vec_numSigmaPoints_;        
    };
}

#endif // MONTE_CARLO_PROBLEM_DESCRIPTION_HPP
