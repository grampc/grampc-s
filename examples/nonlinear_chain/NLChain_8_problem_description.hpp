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

#ifndef NLCHAIN8_HPP
#define NLCHAIN8_HPP

#include "problem_description/problem_description.hpp"

using namespace grampc;

class Chain_8_ProblemDescription : public grampc::ProblemDescription
{
public:
    Chain_8_ProblemDescription(const std::vector<typeRNum>& pCost);

    virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;

    virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;

    virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
    virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
    virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;

    virtual void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
    virtual void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
	virtual void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;

    virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
    virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;

    virtual void hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) override;
    virtual void dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
    virtual void dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;

private:
    // cost parameters
    std::vector<typeRNum> pCost_;
};

#endif // NLCHAIN8_HPP