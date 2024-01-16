#include "reactor_problem_description.hpp"

ReactorProblemDescription::ReactorProblemDescription(const std::vector<typeRNum>& pSys,
                                                     const std::vector<typeRNum>& pCost,
                                                     const std::vector<typeRNum>& pCon)
    : pSys_(pSys), pCost_(pCost), pCon_(pCon)
{
}

void ReactorProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = 2;
    *Nu = 1;
    *Np = 0;
    *Ng = 0;
    *Nh = 1;
    *NgT = 0;
    *NhT = 0;
}

void ReactorProblemDescription::ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = -pSys_[0] * x[0] - pSys_[2] * x[0] * x[0] + (1 - x[0]) * u[0];
    out[1] = pSys_[0] * x[0] - pSys_[1] * x[1] - x[1] * u[0];
}

void ReactorProblemDescription::dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = (-pSys_[0] - pSys_[2] * 2 * x[0] - u[0]) * adj[0] + pSys_[0] * adj[1];
    out[1] = (-pSys_[1] - u[0]) * adj[1];
}

void ReactorProblemDescription::dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *adj, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = (1 - x[0]) * adj[0] + (-x[1]) * adj[1];
}

void ReactorProblemDescription::lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = pCost_[2] * (x[0] - xdes[0]) * (x[0] - xdes[0]) +
             pCost_[3] * (x[1] - xdes[1]) * (x[1] - xdes[1]) +
             pCost_[4] * (u[0] - udes[0]) * (u[0] - udes[0]);
}

void ReactorProblemDescription::dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = 2 * pCost_[2] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[3] * (x[1] - xdes[1]);
}

void ReactorProblemDescription::dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes)
{
    out[0] = 2 * pCost_[4] * (u[0] - udes[0]);
}

void ReactorProblemDescription::Vfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = pCost_[0] * (x[0] - xdes[0]) * (x[0] - xdes[0]) +
             pCost_[1] * (x[1] - xdes[1]) * (x[1] - xdes[1]);
}

void ReactorProblemDescription::dVdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes)
{
    out[0] = 2 * pCost_[0] * (x[0] - xdes[0]);
    out[1] = 2 * pCost_[1] * (x[1] - xdes[1]);
}

void ReactorProblemDescription::hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p)
{
    out[0] = x[1] - pCon_[0];
}

void ReactorProblemDescription::dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0.0;
    out[1] = vec[0];
}

void ReactorProblemDescription::dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec)
{
    out[0] = 0.0;
}
