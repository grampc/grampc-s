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


#ifndef RESAMPLING_GP_PROBLEM_DESCRIPTION_HPP
#define RESAMPLING_GP_PROBLEM_DESCRIPTION_HPP

#include "point_transformation/point_transformation.hpp"
#include "constraint_approx/chebyshev_constraint_approximation.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"
#include "gaussian_process/gaussian_process.hpp"

namespace grampc
{
    // Interface to GRAMPC using a point-based approximation of the uncertainties with a GP-model of unknown parts of the system dynamics
    class ResamplingGPProblemDescription : public ProblemDescription
    {
    public:
        // Constructor specifying the approximation type of the constraints and the type of the point-based uncertainty approximation
        ResamplingGPProblemDescription(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                    PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP);

        // Constructor specifying the type of the point-based uncertainty approximation
        ResamplingGPProblemDescription(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation, 
                                    const std::vector<GaussianProcessPtr>& gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP);        

        /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
			inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
        virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;


        /** System function f(t,x,u,p,userparam)
		------------------------------------ **/
        virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
        
        // Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx)
        virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) override;
        
        // Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du)
        virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) override;


        /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
        virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
        
        // Gradient dl/dx
        virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
        
        // Gradient dl/du
        virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;


        /** Terminal cost V(T,x,p)
		-------------------------------------------------- **/
        virtual void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
        
        // Gradient dV/dx
        virtual void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;

        /** Gradient dV/dT **/
		virtual void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;


        /** Inequality constraints h(t,x,u,p) < 0
		-------------------------------------------------- **/
        virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
        
        //  Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx)
        virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;

        // Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du)
        virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;


        /** Terminal inequality constraints hT(T,x,p) < 0 
         * ------------------------------------------------ **/
        virtual void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) override;

		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
        
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;


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
        const std::vector<GaussianProcessPtr> gaussianProcessVec_;
        const std::vector<typeInt> dynamicsIndicesWithGP_;
        Vector x0_;
        Vector p0_;
        Vector constraintStdDev_;
        Vector terminalConstraintStdDev_;
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
	typedef std::shared_ptr<ResamplingGPProblemDescription> ResamplingGPProblemDescriptionPtr;

    // Constructor specifying the approximation type of the constraints and the type of the point-based uncertainty approximation
	ResamplingGPProblemDescriptionPtr ResamplingGPProblem(ProblemDescriptionPtr problemDescription, ChanceConstraintApproximationConstPtr constraintApproximation,
                                            PointTransformationPtr pointTransformation, const std::vector<GaussianProcessPtr>& gaussianProcessVec,
                                            const std::vector<typeInt>& dynamicsIndicesWithGP);
                                
    // Constructor specifying the type of the point-based uncertainty approximation
    ResamplingGPProblemDescriptionPtr ResamplingGPProblem(ProblemDescriptionPtr problemDescription, PointTransformationPtr pointTransformation, 
                                            const std::vector<GaussianProcessPtr>& gaussianProcessVec, const std::vector<typeInt>& dynamicsIndicesWithGP);
}

#endif // RESAMPLING_GP_PROBLEM_DESCRIPTION_HPP
