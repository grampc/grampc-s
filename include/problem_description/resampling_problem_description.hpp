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


#ifndef RESAMPLING_PROBLEM_DESCRIPTION_HPP
#define RESAMPLING_PROBLEM_DESCRIPTION_HPP

#include "point_transformation/point_transformation.hpp"
#include "constraint_approx/chance_constraint_approximation.hpp"
#include "constraint_approx/chebyshev_constraint_approximation.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"

namespace grampc
{
    // Interface to GRAMPC using a point-based approximation of the uncertainties
    class ResamplingProblemDescription : public ProblemDescription
    {
    public:
        // Constructor specifying the approximation type of the constraints and the type of the point-based uncertainty approximation
        ResamplingProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation);
        
        // Constructor with zero-mean Wiener process and specified diffusion matrix
        ResamplingProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation, MatrixConstRef diffMatrixWienerProcess);

        // Constructor specifying the type of the point-based uncertainty approximation
        ResamplingProblemDescription(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation);
        
        // Constructor with zero-mean Wiener process and specified diffusion matrix
        ResamplingProblemDescription(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation, MatrixConstRef diffMatrixWienerProcess);

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
        typeInt numSigmaPoints_;
        typeInt numStates_;
        typeInt numParams_;
        typeInt numControlInputs_;
        typeInt numConstraints_;
        typeInt numTerminalConstraints_;
        typeInt pointDim_;
                
        ProblemDescriptionPtr problemDescription_;
        PointTransformationPtr pointTransformation_;
        Vector constraintTighteningCoeff_;
        Vector x0_;
        Vector p0_;
        Vector constraintStdDev_;
        Vector terminalConstraintStdDev_;
        Matrix diffMatrixWienerProcess_;
        Matrix pointsTransformed_;
        Matrix covStateAndParam_;
        Vector dcov_dPointsY_;
        Vector meanStateAndParam_;
        Matrix constraintMatrix_;
        Matrix terminalConstraintMatrix_;
        RowVector dmean1D_dPoints_;
        Matrix constraintVec_;
        Matrix terminalConstraintVec_;
        Vector tempScalar_;

        // Standard Cholesky decomposition
        Eigen::LLT<Matrix> llt_;

        // Temporary vectors
        Vector temp_vec_pointDim_numPoints_;
        Vector temp_vec_numInputs_;
    };

    // Alias
	typedef std::shared_ptr<ResamplingProblemDescription> ResamplingProblemDescriptionPtr;
	
    // Constructor specifying the approximation type of the constraints and the type of the point-based uncertainty approximation
    ResamplingProblemDescriptionPtr ResamplingProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation);
    // Constructor with zero-mean Wiener process and specified diffusion matrix
    ResamplingProblemDescriptionPtr ResamplingProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation, MatrixConstRef diffMatrixWienerProcess);
    // Constructor specifying the type of the point-based uncertainty approximation
    ResamplingProblemDescriptionPtr ResamplingProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation);
    // Constructor with zero-mean Wiener process and specified diffusion matrix
    ResamplingProblemDescriptionPtr ResamplingProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation, MatrixConstRef diffMatrixWienerProcess);
}

#endif // RESAMPLING_PROBLEM_DESCRIPTION_HPP
