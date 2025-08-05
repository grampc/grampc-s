
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

#ifndef INVERTED_PENDULUM_PROBLEM_DESCRIPTION_HPP
#define INVERTED_PENDULUM_PROBLEM_DESCRIPTION_HPP

#include "problem_description/problem_description.hpp"

using namespace grampc;

#define POW2(x) ((x) * (x))

class InvertedPendulumProblemDescription : public grampc::ProblemDescription
{
public:
    InvertedPendulumProblemDescription(const std::vector<typeRNum> &pSys, const std::vector<typeRNum> &pCost, const std::vector<typeRNum> &pCon);

    // Specification of the dimensions
    virtual void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT) override;

    // System dynamics and derivatives
    virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef adj, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;

    // Integral cost and its derivatives
    virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
    virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
    virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
    virtual void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;

    // Terminal cost and its derivatives
    virtual void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
    virtual void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
    virtual void dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
    virtual void dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;

    // Constraints and its derivatives
    virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
    virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;

    // Additional functions required for Taylor-SMPC
    virtual void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dfdpdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
    virtual void dhdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;

private:
    // system parameters
    std::vector<typeRNum> pSys_;

    // cost parameters
    std::vector<typeRNum> pCost_;

    // constraint parameters
    std::vector<typeRNum> pCon_;

    typeRNum temp_;
};

#endif // INVERTED_PENDULUM_PROBLEM_DESCRIPTION_HPP
