/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */

// define macros
#define S_FUNCTION_NAME  grampc_s_Sfct
#define S_FUNCTION_LEVEL 2

// C headers
#ifdef __cplusplus
extern "C" {
#endif

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#include "simstruc.h"

#ifdef __cplusplus
} // end of extern "C" scope
#endif

// C++ headers
#include "grampc_s.hpp"
#include "grampc_interface/grampc_interface.hpp"
#include "inverted_pendulum_problem_description.hpp"


static void mdlInitializeSizes(SimStruct *S)
{
    typeInt numStates = 4;
    typeInt numInputs = 1;
    typeInt Nhor = 30;

    // index variable
	typeInt  i;

	// s-function parameters
	ssSetNumSFcnParams(S, 0);
	/* parameter tuning not allowed during simulation */
	for (i = 0; i < ssGetSFcnParamsCount(S); i++) {
		ssSetSFcnParamTunable(S, i, SS_PRM_NOT_TUNABLE);
	}

	// sizes of input ports
	if (!ssSetNumInputPorts(S, 3)) return;                     /* number of input ports */
	ssSetInputPortWidth(S, 0, numStates);                       /* no. of inputs at input port 0 (x)       */
	ssSetInputPortWidth(S, 1, 1);                               /* no. of inputs at input port 1 (t)       */
    ssSetInputPortWidth(S, 2, MAX(numStates * numStates, 1));   /* no. of inputs at input port 2 (xCovar)  */

    // set options for input ports
	for (i = 0; i < ssGetNumInputPorts(S); i++) {
		ssSetInputPortDirectFeedThrough(S, i, 1);               /* direct input feedthrough */
		ssSetInputPortRequiredContiguous(S, i, 1);              /* must be contiguous to use  ssGetInputPortRealSignal*/
		ssSetInputPortDataType(S, i, SS_TYPERNUM);
	}

	// sizes of output ports
	if (!ssSetNumOutputPorts(S, 4)) return;                     /* number of output ports                     */
	ssSetOutputPortWidth(S, 0, MAX(numInputs, 1));              /* no. of outputs at output port 0 (unext)    */
    ssSetOutputPortWidth(S, 1, Nhor);                           /* no. of outputs at output port 1 (uPred)    */
    ssSetOutputPortWidth(S, 2, MAX(numInputs, 1));              /* no. of outputs at output port 2 (compTime) */
    ssSetOutputPortWidth(S, 3, 2);                              /* no. of outputs at output port 3 (cost) */

    // set options for output ports
	for (i = 0; i < ssGetNumOutputPorts(S); i++) {
		ssSetOutputPortDataType(S, i, SS_TYPERNUM);
	}

	// other sizes
	ssSetNumSampleTimes(S, 0);    // number of sample times
	ssSetNumRWork(S, 0);          // number of real work vectors
	ssSetNumIWork(S, 0);          // number of integer work vectors
	ssSetNumPWork(S, 2);         // number of pointer work vectors
	ssSetNumModes(S, 0);          // number of modes to switch between
	ssSetNumNonsampledZCs(S, 0);  // number of non-sampled zero-crossings
}



#if defined(MATLAB_MEX_FILE)
#define MDL_SET_INPUT_PORT_WIDTH
//Set the width of an input port that accepts 1-D (vector) signals
static void mdlSetInputPortWidth(SimStruct *S, typeInt port, typeInt inputPortWidth)
{
	if ((ssGetInputPortWidth(S, port) != DYNAMICALLY_SIZED) && (ssGetInputPortWidth(S, port) != inputPortWidth))
	{
		ssSetErrorStatus(S, "something wrong");
		return;
	}
	ssSetInputPortWidth(S, port, inputPortWidth);
}

#define MDL_SET_OUTPUT_PORT_WIDTH
static void mdlSetOutputPortWidth(SimStruct *S, typeInt port, typeInt outputPortWidth)
{
	/* We have no dynamic outputs */
}

#define MDL_SET_DEFAULT_PORT_DIMENSION_INFO
static void mdlSetDefaultPortDimensionInfo(SimStruct *S)
{
	ssSetInputPortWidth(S, 4, 1); /* default to one */
}

#define MDL_SET_WORK_WIDTHS
static void mdlSetWorkWidths(SimStruct *S)
{
	/* We have no states */
}
#endif


/***************************************************************************
 * mdlInitializeSampleTimes                                                *
 ***************************************************************************/

//Specify the sample rates at which this C MEX S-function operates
static void mdlInitializeSampleTimes(SimStruct *S)
{
	/* Attention: Read dt */
	ssSetSampleTime(S, 0, -1);
	ssSetOffsetTime(S, 0, 0.0);
}


/***************************************************************************
 * mdlStart                                                                *
 ***************************************************************************/

#define MDL_START
#if defined (MDL_START)
static void mdlStart(SimStruct *S)
{
    constexpr typeRNum pi = 3.14159265358979323846;

    // dimensions
    typeInt numStates = 4;
    typeInt numInputs = 1;
    typeInt Nhor = 30;
    typeRNum UT_alpha = 0.1;
    typeRNum UT_beta  = 2.0;
    typeRNum UT_kappa = 0;

    grampc::PointTransformationPtr transform = grampc::UT(numStates, numStates, UT_alpha, UT_beta, UT_kappa);

    // parameters
    typeRNum Thor = 0.7;
    typeRNum dt = 0.001;
   
    std::vector<typeRNum> u0 = {0.0};
    std::vector<typeRNum> constraintsAbsTol = {1e-6, 1e-6};

    // desired values for states and inputs
    std::vector<typeRNum> xdes = {0, 0, pi, 0};
    std::vector<typeRNum> udes = {0.0};

    // input constraints
    std::vector<typeRNum> umax = {5};
    std::vector<typeRNum> umin = {-5};

    // declare vectors for the problem description
    std::vector<typeRNum> pSys = {2e-3, 0.2, 0.15, 9.81, 1.3e-3};
    std::vector<typeRNum> pCost = {100, 1, 1, 1, 0, 0, 0, 0, 1e-9};
    std::vector<typeRNum> pCon = {0.8};

    // Chance constratypeInt probability
    grampc::Vector vec_chance_constraint(2);
    vec_chance_constraint << 0.95, 0.95, 0.95;

    // Chance constratypeInt approximation
    grampc::ChanceConstraintApproximationPtr constraintApprox = grampc::GaussianApprox(vec_chance_constraint);

    // create nominal problem description
    grampc::ProblemDescriptionPtr pendulumProblem = grampc::ProblemDescriptionPtr(new InvertedPendulumProblemDescription(pSys, pCost, pCon));

    // configure stochastic problem description
    ssGetPWork(S)[0] = (void *) new grampc::SigmaPointProblemDescriptionPtr(new grampc::SigmaPointProblemDescription(pendulumProblem, constraintApprox, transform));
    grampc::SigmaPointProblemDescriptionPtr *problem = (grampc::SigmaPointProblemDescriptionPtr *) ssGetPWork(S)[0];

    // create solver
    ssGetPWork(S)[1] = (void *) new grampc::Grampc(*problem);
    grampc::Grampc *solver = (grampc::Grampc *) ssGetPWork(S)[1];

    const typeGRAMPCparam* par = solver->getParameters();
    xdes.resize(par->Nx, 0.0);
    constraintsAbsTol.resize(par->Nc, 0.0);

    // set parameters
    solver->setparam_real("Thor", Thor);
    solver->setparam_real("dt", dt);
    solver->setparam_real_vector("xdes", &xdes[0]);
    solver->setparam_real_vector("u0", &u0[0]);
    solver->setparam_real_vector("udes", &udes[0]);
    solver->setparam_real_vector("umax", &umax[0]);
    solver->setparam_real_vector("umin", &umin[0]);

    // set options
    solver->setopt_int("MaxGradIter", 3);
    solver->setopt_int("MaxMultIter", 2);
    solver->setopt_int("Nhor", 20);
    solver->setopt_string("Integrator", "heun");
    solver->setopt_real("PenaltyMin", 1e3);
    solver->setopt_real_vector("ConstraintsAbsTol", &constraintsAbsTol[0]);
}
#endif   /*MDL_START*/


/***************************************************************************
 * mdlOutputs                                                              *
 ***************************************************************************/

static void mdlOutputs(SimStruct *S, typeInt tid)
{
    std::chrono::steady_clock::time_point ticFunctionStart = std::chrono::steady_clock::now();

    // index variable
	typeInt  i;

    // get pointers to objects defined in mdlStart()
    grampc::SigmaPointProblemDescriptionPtr *problem = (grampc::SigmaPointProblemDescriptionPtr *) ssGetPWork(S)[0];
    grampc::Grampc *grampc = (grampc::Grampc *) ssGetPWork(S)[1];

	// pointers to the inputs of the s-function
	ctypeRNum *x0 = (typeRNum *)ssGetInputPortRealSignal(S, 0);
	ctypeRNum *t0 = (typeRNum *)ssGetInputPortRealSignal(S, 1);
    ctypeRNum *xCovar = (typeRNum *)ssGetInputPortRealSignal(S, 2);

    // pointers to outputs of the s-function
	typeRNum *unext = (typeRNum *)ssGetOutputPortRealSignal(S, 0);
    typeRNum *uPred = (typeRNum *)ssGetOutputPortRealSignal(S, 1);
    typeRNum *compTime = (typeRNum *)ssGetOutputPortRealSignal(S, 2);
    typeRNum *cost = (typeRNum *)ssGetOutputPortRealSignal(S, 3);
  
    // dimensions
    typeInt numStates = 4;
    typeInt numInputs = 1;
    typeInt Nhor = 30;

    // set mean and covariance of states
    grampc::DistributionPtr state = grampc::Gaussian(Eigen::Map<const grampc::Vector>(x0,numStates), Eigen::Map<const grampc::Matrix>(xCovar, numStates, numStates));

    // remove warning unreferenced formal parameter
	(void)(tid);

	// set time, states, and parameters based on s-function inputs
	grampc->setparam_real("t0", t0[0]);
    (*problem)->compute_x0_and_p0(state);
    grampc->setparam_real_vector("x0", (*problem)->x0());
    grampc->setparam_real_vector("p0", (*problem)->p0());

	// run solver
	grampc->run();
    std::chrono::steady_clock::time_point toc = std::chrono::steady_clock::now();

    // Computation time in millisecond
    std::chrono::duration<double, std::milli> timeTotal_ms = toc - ticFunctionStart;

    // solution structure and workspace structure
	const typeGRAMPCsol *sol = grampc->getSolution();
    const typeGRAMPCrws *rws = grampc->getWorkspace();

	// Copy the solution to the outputs of the s-function
	MatCopy(unext, sol->unext, 1, numInputs);
    MatCopy(uPred, rws->u, 1, Nhor);
    compTime[0] = timeTotal_ms.count();
    MatCopy(cost, sol->J, 1, 2); 
}


/****************************************************************************
 * mdlTerminate                                                             *
 ****************************************************************************/
static void mdlTerminate(SimStruct *S)
{
    // get pointers to objects defined in mdlStart()
    grampc::SigmaPointProblemDescriptionPtr *problem = (grampc::SigmaPointProblemDescriptionPtr *) ssGetPWork(S)[0];
    grampc::Grampc *grampc = (grampc::Grampc *) ssGetPWork(S)[1];

    // free allocated memory
    delete problem;
	delete grampc;
}


/****************************************************************************
 * Compiler dependent settings                                              *
 ****************************************************************************/

#ifdef	MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif

 /****************************************************************************
	* End of C-Code S-Function	                                               *
	****************************************************************************/
