
#include "problem_description/grampc_interface.hpp"


namespace grampc
{
	Grampc::Grampc(StochasticProblemDescription *problem_description)
		: problem_description_(problem_description)
	{
#ifdef FIXEDSIZE
        grampc_struct_.param = &grampc_param_;
        grampc_struct_.opt = &grampc_opt_;
        grampc_struct_.sol = &grampc_sol_;
        grampc_struct.rws = &grampc_rws_;
        grampc_ = &grampc_struct;
#endif
		grampc_init(&grampc_, problem_description_);
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

	void Grampc::setopt_real_vector(const char *optName, ctypeRNum *optValue)
	{
		grampc_setopt_real_vector(grampc_, optName, optValue);
	}

	void Grampc::setopt_int_vector(const char *optName, ctypeInt *optValue)
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

	void Grampc::setparam_real_vector(const char *paramName, ctypeRNum *paramValue)
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

}
