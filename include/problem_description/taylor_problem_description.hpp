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


#ifndef TAYLOR_PROBLEM_DESCRIPTION_HPP
#define TAYLOR_PROBLEM_DESCRIPTION_HPP

#include "point_transformation/point_transformation.hpp"
#include "constraint_approx/chance_constraint_approximation.hpp"
#include "constraint_approx/chebyshev_constraint_approximation.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    // Interface to GRAMPC using a first order Taylor series approximation for uncertainty propagation
    class TaylorProblemDescription : public ProblemDescription
    {
    public:
        // Constructor specifying the type of constraint approximation
        TaylorProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation);

        // Constructor with zero-mean Wiener process and specified diffusion matrix
        TaylorProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef diffMatrixWienerProcess);

        // Constructor for problems without constraints
        TaylorProblemDescription(ProblemDescriptionPtr problemDescription);

        // Constructor with zero-mean Wiener process and specified diffusion matrix
        TaylorProblemDescription(ProblemDescriptionPtr problemDescription, MatrixConstRef diffMatrixWienerProcess);
        
        /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
			inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;


        /** System function f(t,x,u,p,userparam)
		------------------------------------ **/
        virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param) override;
        
        // Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx)
        virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param) override;
        
        // Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du)
        virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef adj, const typeGRAMPCparam *param) override;


        /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
        virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param) override;
        
        // Gradient dl/dx
        virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param) override;
        
        // Gradient dl/du
        virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param) override;


        /** Terminal cost V(T,x,p)
		-------------------------------------------------- **/
        virtual void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param) override;
        
        // Gradient dV/dx
        virtual void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param) override;

        /** Gradient dV/dT **/
		virtual void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param) override;


        /** Inequality constraints h(t,x,u,p) < 0
		-------------------------------------------------- **/
        virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param) override;
        
        //  Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx)
        virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param) override;

        // Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du)
        virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param) override;


        /** Terminal inequality constraints hT(T,x,p) < 0 
         * ------------------------------------------------ **/
        virtual void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param) override;

		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param) override;
        
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param) override;


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
        
        ProblemDescriptionPtr problemDescription_;
        Vector constraintTighteningCoeff_;
        Matrix diffMatrixWienerProcess_;
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

    // Alias
	typedef std::shared_ptr<TaylorProblemDescription> TaylorProblemDescriptionPtr;
	
    // Constructor specifying the type of constraint approximation
    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation);
    // Constructor with zero-mean Wiener process and specified diffusion matrix
    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation, MatrixConstRef diffMatrixWienerProcess);
    // Constructor for problems without constraints
    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription);
    // Constructor with zero-mean Wiener process and specified diffusion matrix
    TaylorProblemDescriptionPtr TaylorProblem(ProblemDescriptionPtr problemDescription, MatrixConstRef diffMatrixWienerProcess);
}

#endif // TAYLOR_PROBLEM_DESCRIPTION_HPP
