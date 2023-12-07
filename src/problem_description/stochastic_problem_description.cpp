#include "problem_description/stochastic_problem_description.hpp"


extern "C"
{
	void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->ocp_dim(Nx, Nu, Np, Ng, Nh, NgT, NhT);
	}


	void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->ffct(out, t, x, u, p);
	}
	void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dfdx_vec(out, t, x, adj, u, p);
	}
	void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dfdu_vec(out, t, x, adj, u, p);
	}
	void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dfdp_vec(out, t, x, adj, u, p);
	}


	void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->lfct(out, t, x, u, p, xdes, udes);
	}
	void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dldx(out, t, x, u, p, xdes, udes);
	}
	void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dldu(out, t, x, u, p, xdes, udes);
	}
	void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dldp(out, t, x, u, p, xdes, udes);
	}


	void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->Vfct(out, t, x, p, xdes);
	}
	void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dVdx(out, t, x, p, xdes);
	}
	void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dVdp(out, t, x, p, xdes);
	}
	void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dVdT(out, t, x, p, xdes);
	}


	void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->gfct(out, t, x, u, p);
	}
	void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgdx_vec(out, t, x, u, p, vec);
	}
	void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgdu_vec(out, t, x, u, p, vec);
	}
	void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgdp_vec(out, t, x, u, p, vec);
	}


	void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->hfct(out, t, x, u, p);
	}
	void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhdx_vec(out, t, x, u, p, vec);
	}
	void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhdu_vec(out, t, x, u, p, vec);
	}
	void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhdp_vec(out, t, x, u, p, vec);
	}


	void gTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->gTfct(out, t, x, p);
	}
	void dgTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgTdx_vec(out, t, x, p, vec);
	}
	void dgTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgTdp_vec(out, t, x, p, vec);
	}
	void dgTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dgTdT_vec(out, t, x, p, vec);
	}


	void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->hTfct(out, t, x, p);
	}
	void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdx_vec(out, t, x, p, vec);
	}
	void dhTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdp_vec(out, t, x, p, vec);
	}
	void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdT_vec(out, t, x, p, vec);
	}


	void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dhTdT_vec(out, t, x, u, p);
	}
	void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *adj, ctypeRNum *p, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->dHdxdt(out, t, x, u, adj, p);
	}
	void Mfct(typeRNum *out, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->Mfct(out);
	}
	void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
	{
		((grampc::StochasticProblemDescription*)userparam)->Mtrans(out);
	}

}
