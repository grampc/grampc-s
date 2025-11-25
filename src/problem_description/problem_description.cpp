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


#include "problem_description/problem_description.hpp"

extern "C"
{
	void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
	{
		((grampc::ProblemDescription*)userparam)->ocp_dim(Nx, Nu, Np, Ng, Nh, NgT, NhT);
	}

	void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->ffct(outMap, t, xMap, uMap, pMap, param);
	}

	void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *adj, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> adjMap(adj, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->dfdx_vec(outMap, t, xMap, uMap, pMap, adjMap, param);
	}

	void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *adj, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> adjMap(adj, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nu_);

		((grampc::ProblemDescription*)userparam)->dfdu_vec(outMap, t, xMap, uMap, pMap, adjMap, param);
	}

	void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *adj, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> adjMap(adj, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Np_);

		((grampc::ProblemDescription*)userparam)->dfdp_vec(outMap, t, xMap, uMap, pMap, adjMap, param);
	}


	void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, 1);

		((grampc::ProblemDescription*)userparam)->lfct(outMap, t, xMap, uMap, pMap, param);
	}

	void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->dldx(outMap, t, xMap, uMap, pMap, param);
	}

	void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nu_);

		((grampc::ProblemDescription*)userparam)->dldu(outMap, t, xMap, uMap, pMap, param);
	}

	void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Np_);

		((grampc::ProblemDescription*)userparam)->dldp(outMap, t, xMap, uMap, pMap, param);
	}


	void Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, 1);

		((grampc::ProblemDescription*)userparam)->Vfct(outMap, t, xMap, pMap, param);
	}

	void dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->dVdx(outMap, t, xMap, pMap, param);
	}

	void dVdp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Np_);

		((grampc::ProblemDescription*)userparam)->dVdp(outMap, t, xMap, pMap, param);
	}

	void dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, 1);

		((grampc::ProblemDescription*)userparam)->dVdT(outMap, t, xMap, pMap, param);
	}


	void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}

	void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}

	void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}


	void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nh_);

		((grampc::ProblemDescription*)userparam)->hfct(outMap, t, xMap, uMap, pMap, param);
	}

	void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->Nh_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->dhdx_vec(outMap, t, xMap, uMap, pMap, vecMap, param);
	}

	void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->Nh_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nu_);

		((grampc::ProblemDescription*)userparam)->dhdu_vec(outMap, t, xMap, uMap, pMap, vecMap, param);
	}

	void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> uMap(u, ((grampc::ProblemDescription*)userparam)->Nu_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->Nh_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Np_);

		((grampc::ProblemDescription*)userparam)->dhdp_vec(outMap, t, xMap, uMap, pMap, vecMap, param);
	}


	void gTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}

	void dgTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}

	void dgTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}

	void dgTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
	}


	void hTfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->NhT_);

		((grampc::ProblemDescription*)userparam)->hTfct(outMap, t, xMap, pMap, param);
	}

	void dhTdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->NhT_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Nx_);

		((grampc::ProblemDescription*)userparam)->dhTdx_vec(outMap, t, xMap, pMap, vecMap, param);
	}

	void dhTdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->NhT_);
		Eigen::Map<grampc::Vector> outMap(out, ((grampc::ProblemDescription*)userparam)->Np_);

		((grampc::ProblemDescription*)userparam)->dhTdp_vec(outMap, t, xMap, pMap, vecMap, param);
	}

	void dhTdT_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, const typeGRAMPCparam *param, typeUSERPARAM* userparam)
	{
		Eigen::Map<const grampc::Vector> xMap(x, ((grampc::ProblemDescription*)userparam)->Nx_);
		Eigen::Map<const grampc::Vector> pMap(p, ((grampc::ProblemDescription*)userparam)->Np_);
		Eigen::Map<const grampc::Vector> vecMap(vec, ((grampc::ProblemDescription*)userparam)->NhT_);
		Eigen::Map<grampc::Vector> outMap(out, 1);

		((grampc::ProblemDescription*)userparam)->dhTdT_vec(outMap, t, xMap, pMap, vecMap, param);
	}

	void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *adj, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void Mfct(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}

	void Mtrans(typeRNum *out, const typeGRAMPCparam *param, typeUSERPARAM *userparam)
	{
	}
}
