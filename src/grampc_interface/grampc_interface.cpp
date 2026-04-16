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


#include "grampc_interface/grampc_interface.hpp"
#include "python/pygrampc_types.hpp"

namespace grampc
{

    GrampcParam::GrampcParam()
        : x0(NULL, 0), xdes(NULL, 0), u0(NULL, 0), udes(NULL, 0),
        umax(NULL, 0), umin(NULL, 0), p0(NULL, 0), pmax(NULL, 0), pmin(NULL, 0)
    {}

    GrampcParam::GrampcParam(const typeGRAMPCparam* param)
        : Nx(&param->Nx),
          Nu(&param->Nu),
          Np(&param->Np),
          Ng(&param->Ng),
          Nh(&param->Nh),
          NgT(&param->NgT),
          NhT(&param->NhT),
          Nc(&param->Nc),
          x0(param->x0, *Nx), 
          xdes(param->xdes, *Nx), 
          u0(param->u0, *Nu), 
          udes(param->udes, *Nu),
          umax(param->umax, *Nu), 
          umin(param->umin, *Nu), 
          p0(param->p0, *Np), 
          pmax(param->pmax, *Np), 
          pmin(param->pmin, *Np),
          Thor(&param->Thor),
          Tmax(&param->Tmax),
          Tmin(&param->Tmin),
          dt(&param->dt),
          t0(&param->t0)
    {}

    void GrampcParam::remap_memory(const typeGRAMPC *grampc)
    {
        Nx = &grampc->param->Nx;
        Nu = &grampc->param->Nu;
        Np = &grampc->param->Np;
        Ng = &grampc->param->Ng;
        Nh = &grampc->param->Nh;
        NgT = &grampc->param->NgT;
        NhT = &grampc->param->NhT;
        Nc = &grampc->param->Nc;

        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&x0) Eigen::Map<Vector> (grampc->param->x0, grampc->param->Nx);
        new (&xdes) Eigen::Map<Vector>(grampc->param->xdes, grampc->param->Nx);

        new (&u0) Eigen::Map<Vector>(grampc->param->u0, grampc->param->Nu);
        new (&udes) Eigen::Map<Vector>(grampc->param->udes, grampc->param->Nu);
        new (&umax) Eigen::Map<Vector>(grampc->param->umax, grampc->param->Nu);
        new (&umin) Eigen::Map<Vector>(grampc->param->umin, grampc->param->Nu);

        new (&p0) Eigen::Map<Vector>(grampc->param->p0, grampc->param->Np);
        new (&pmax) Eigen::Map<Vector>(grampc->param->pmax, grampc->param->Np);
        new (&pmin) Eigen::Map<Vector>(grampc->param->pmin, grampc->param->Np);

        Thor = &grampc->param->Thor;
        Tmax = &grampc->param->Tmax;
        Tmin = &grampc->param->Tmin;

        dt = &grampc->param->dt;
        t0 = &grampc->param->t0;
    }

	Grampc::Grampc(ProblemDescriptionPtr problem_description)
		: problem_description_(problem_description)
	{
#ifdef FIXEDSIZE
        grampc_struct_.param = &grampc_param_;
        grampc_struct_.opt = &grampc_opt_;
        grampc_struct_.sol = &grampc_sol_;
        grampc_struct.rws = &grampc_rws_;
        grampc_ = &grampc_struct;
#endif
		grampc_init(&grampc_, problem_description_.get());		
	}

	Grampc::~Grampc()
	{
		grampc_free(&grampc_);
	}


	void Grampc::setopt_real(const char *optName, ctypeRNum optValue)
	{
		grampc_setopt_real(grampc_, optName, optValue);
	}

	void Grampc::setopt_int(const char *optName, ctypeInt optValue)
	{
		grampc_setopt_int(grampc_, optName, optValue);
	}

	void Grampc::setopt_string(const char *optName, const char *optValue)
	{
		grampc_setopt_string(grampc_, optName, optValue);
	}

	void Grampc::setopt_real_vector(const char *optName, VectorConstRef optValue)
	{
		grampc_setopt_real_vector(grampc_, optName, optValue.data());
	}
	void Grampc::setopt_real_vector(const char *optName, ctypeRNum* optValue)
	{
		grampc_setopt_real_vector(grampc_, optName, optValue);
	}

	void Grampc::setopt_int_vector(const char *optName, IntVectorConstRef optValue)
	{
		grampc_setopt_int_vector(grampc_, optName, optValue.data());
	}
	void Grampc::setopt_int_vector(const char *optName, ctypeInt* optValue)
	{
		grampc_setopt_int_vector(grampc_, optName, optValue);
	}

	void Grampc::printopt() const
	{
		grampc_printopt(grampc_);
	}


	void Grampc::setparam_real(const char *paramName, ctypeRNum paramValue)
	{
		grampc_setparam_real(grampc_, paramName, paramValue);
	}

	void Grampc::setparam_real_vector(const char *paramName, VectorConstRef paramValue)
	{
		grampc_setparam_real_vector(grampc_, paramName, paramValue.data());
	}
	void Grampc::setparam_real_vector(const char *paramName, ctypeRNum* paramValue)
	{
		grampc_setparam_real_vector(grampc_, paramName, paramValue);
	}

	void Grampc::printparam() const
	{
		grampc_printparam(grampc_);
	}


    typeInt Grampc::estim_penmin(ctypeInt rungrampc)
    {
        return grampc_estim_penmin(grampc_, rungrampc);
    }


	void Grampc::run()
	{
		grampc_run(grampc_);
	}

	typeInt Grampc::printstatus(ctypeInt status, ctypeInt level)
	{
		return grampc_printstatus(status, level);
	}

	const typeGRAMPCopt* Grampc::getOptions() const
	{
		return grampc_->opt;
	}

	const typeGRAMPCparam* Grampc::getParameters() const
	{
		return grampc_->param;
	}

	const typeGRAMPCrws* Grampc::getWorkspace() const
	{
		return grampc_->rws;
	}

	const typeGRAMPCsol* Grampc::getSolution() const
	{
		return grampc_->sol;
	}

	ProblemDescriptionPtr Grampc::getProblem() const
	{
		return problem_description_;
	}

	GrampcPtr Solver(ProblemDescriptionPtr problem_description)
	{
		return GrampcPtr(new Grampc(problem_description));
	}
}
