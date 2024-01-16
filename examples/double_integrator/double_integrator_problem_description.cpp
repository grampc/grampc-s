#include "double_integrator_problem_description.hpp"

DoubleIntegratorProblemDescription::DoubleIntegratorProblemDescription(const std::vector<typeRNum>& pCost, const std::vector<typeRNum>& pCon)
    : pCost_(pCost), pCon_(pCon)
{
}

void DoubleIntegratorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 2;
	*Nu = 1;
	*Np = 0;
	*Nh = 1;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void DoubleIntegratorProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = x[1];
	out[1] = u[0];
}

void DoubleIntegratorProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = 0;
	out[1] = vec[0];
}

void DoubleIntegratorProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = vec[1];
}

void DoubleIntegratorProblemDescription::dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p)
{
}

void DoubleIntegratorProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = pCost_[0] * (u[0] - udes[0]) * (u[0] - udes[0]) + pCost_[1] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[2] * (x[1] - xdes[1]) * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = 2 * pCost_[1] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[2] * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
	out[0] = 2 * pCost_[0] * (u[0] - udes[0]);
}

void DoubleIntegratorProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = pCost_[3] * (x[0] - xdes[0]) * (x[0] - xdes[0]) + pCost_[4] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +  pCost_[5] * t;
}

void DoubleIntegratorProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = 2 * pCost_[3] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[4] * (x[1] - xdes[1]);
}

void DoubleIntegratorProblemDescription::dVdT(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = pCost_[5];
}

void DoubleIntegratorProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = -x[1] + pCon_[0]; 
}

void DoubleIntegratorProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0;
    out[1] = -vec[0];
}

void DoubleIntegratorProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0;
}

void DoubleIntegratorProblemDescription::dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
}