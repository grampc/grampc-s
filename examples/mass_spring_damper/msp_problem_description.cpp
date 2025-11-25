/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 * 
 * This probfct-file describes a modified mass-spring-damper problem from
 * Kapernick, B.: Gradient-Based Nonlinear Model Predictive Control With
 * Constraint Transformation for Fast Dynamical Systems. Dissertation,
 * Ulm University. Shaker, Aachen, Germany (2016)
 * 
 * Based on GRAMPC:
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * 
 */

#include "msp_problem_description.hpp"

/* square macro */
#define POW2(a) ((a)*(a))

MassSpringDamperProblemDescription::MassSpringDamperProblemDescription(typeInt numberOfMasses, const std::vector<typeRNum>& pSys, const std::vector<typeRNum>& pCost)
: numberOfMasses_(numberOfMasses),
  Nx_(2 * numberOfMasses),
  pSys_(pSys),
  pCost_(pCost)
{

}

void MassSpringDamperProblemDescription::ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT)
{
    *Nx = Nx_;
	*Nu = 2;
	*Np = 2;
	*Nh = 0;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

void MassSpringDamperProblemDescription::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
	for (typeInt k = 0; k <= numberOfMasses_ - 1; k++)
	{
		out[k] = x[numberOfMasses_ + k];
	}

	out[numberOfMasses_] = -2 * p[0] / pSys_[0] * x[0] + p[0] / pSys_[0] * x[1] - 2 * p[1] / pSys_[0] * x[numberOfMasses_] + p[1] / pSys_[0] * x[numberOfMasses_ + 1] + 1 / pSys_[0] * u[0];

	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[numberOfMasses_ + k] = p[0] / pSys_[0] * x[k - 1] - 2 * p[0] / pSys_[0] * x[k] + p[0] / pSys_[0] * x[k + 1] + p[1] / pSys_[0] * x[numberOfMasses_ + k - 1] - 2 * p[1] / pSys_[0] * x[numberOfMasses_ + k] + p[1] / pSys_[0] * x[numberOfMasses_ + k + 1];
	}

	out[2 * numberOfMasses_ - 1] = p[0] / pSys_[0] * x[numberOfMasses_ - 2] - 2 * p[0] / pSys_[0] * x[numberOfMasses_ - 1] + p[1] / pSys_[0] * x[2 * numberOfMasses_ - 2] - 2 * p[1] / pSys_[0] * x[2 * numberOfMasses_ - 1] - 1 / pSys_[0] * u[1];
}

void MassSpringDamperProblemDescription::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
	out[0] = -2 * p[0] / pSys_[0] * vec[numberOfMasses_] + p[0] / pSys_[0] * vec[numberOfMasses_ + 1];

	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[k] = p[0] / pSys_[0] * vec[numberOfMasses_ + k - 1] - 2 * p[0] / pSys_[0] * vec[numberOfMasses_ + k] + p[0] / pSys_[0] * vec[numberOfMasses_ + k + 1];
	}

	out[numberOfMasses_ - 1] = p[0] / pSys_[0] * vec[2 * numberOfMasses_ - 2] - 2 * p[0] / pSys_[0] * vec[2 * numberOfMasses_ - 1];

	out[numberOfMasses_] = vec[0] - 2 * p[1] / pSys_[0] * vec[numberOfMasses_] + p[1] / pSys_[0] * vec[numberOfMasses_ + 1];

	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[numberOfMasses_ + k] = vec[k] + p[1] / pSys_[0] * vec[numberOfMasses_ + k - 1] - 2 * p[1] / pSys_[0] * vec[numberOfMasses_ + k] + p[1] / pSys_[0] * vec[numberOfMasses_ + k + 1];
	}

	out[2 * numberOfMasses_ - 1] = vec[numberOfMasses_ - 1] + p[1] / pSys_[0] * vec[2 * numberOfMasses_ - 2] - 2 * p[1] / pSys_[0] * vec[2 * numberOfMasses_ - 1];
}

void MassSpringDamperProblemDescription::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
	out[0] = 1 / pSys_[0] * vec[numberOfMasses_];
	out[1] = -1 / pSys_[0] * vec[2 * numberOfMasses_ - 1];
}

void MassSpringDamperProblemDescription::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
	out[0] = -2.0 / pSys_[0] * x[0] + 1.0 / pSys_[0] * x[1] * vec[numberOfMasses_];
	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[0] += vec[numberOfMasses_ + k] / pSys_[0] * x[k - 1] - 2 / pSys_[0] * x[k] + 1.0 / pSys_[0] * x[k + 1];
	}
	out[0] += vec[2 * numberOfMasses_ - 1] / pSys_[0] * x[numberOfMasses_ - 2] - 2 / pSys_[0] * x[numberOfMasses_ - 1];
	
	
	out[1] = - 2 / pSys_[0] * x[numberOfMasses_] + 1.0 / pSys_[0] * x[numberOfMasses_ + 1] * vec[numberOfMasses_];
	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[1] += vec[numberOfMasses_ + k] / pSys_[0] * x[numberOfMasses_ + k - 1] - 2 / pSys_[0] * x[numberOfMasses_ + k] + 1.0 / pSys_[0] * x[numberOfMasses_ + k + 1];
	}
	out[1] += vec[2 * numberOfMasses_ - 1] / pSys_[0] * x[2 * numberOfMasses_ - 2] - 2 / pSys_[0] * x[2 * numberOfMasses_ - 1];
}

void MassSpringDamperProblemDescription::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
    ctypeRNum *udes = param->udes;
	out[0] = 0;
	for (typeInt i = 0; i <= Nx_ - 1; i++)
	{
		out[0] += pCost_[i] * POW2(x[i] - xdes[i]);
	}
	for (typeInt i = 0; i <= 2 - 1; i++)
	{
		out[0] += pCost_[Nx_ + i] * POW2(u[i] - udes[i]);
	}
	out[0] *= 0.5;
}

void MassSpringDamperProblemDescription::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
	for (typeInt i = 0; i <= Nx_ - 1; i++)
	{
		out[i] = pCost_[i] * (x[i] - xdes[i]);
	}
}

void MassSpringDamperProblemDescription::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *udes = param->udes;
	for (typeInt i = 0; i <= 2 - 1; i++)
	{
		out[i] = pCost_[Nx_ + i] * (u[i] - udes[i]);
	}
}

void MassSpringDamperProblemDescription::Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
	out[0] = 0;
	for (typeInt i = 0; i <= Nx_ - 1; i++)
	{
		out[0] += pCost_[Nx_ + 2 + i] * POW2(x[i] - xdes[i]);
	}
	out[0] *= 0.5;
}

void MassSpringDamperProblemDescription::dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
    ctypeRNum *xdes = param->xdes;
	for (typeInt i = 0; i <= Nx_ - 1; i++)
	{
		out[i] = pCost_[Nx_ + 2 + i] * (x[i] - xdes[i]);
	}
}

void MassSpringDamperProblemDescription::dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, const typeGRAMPCparam *param)
{
}

void MassSpringDamperProblemDescription::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
}

void MassSpringDamperProblemDescription::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}

void MassSpringDamperProblemDescription::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}

void MassSpringDamperProblemDescription::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec, const typeGRAMPCparam *param)
{
}

void MassSpringDamperProblemDescription::dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, const typeGRAMPCparam *param)
{
	for(typeInt i = 0; i < numberOfMasses_; ++i)
	{
		out[i] = 0.0;
	}

	out[numberOfMasses_] = -2 * p[0] / pSys_[0];
	out[numberOfMasses_ + 1] = p[0] / pSys_[0]; 

	for(typeInt i = numberOfMasses_ + 2; i < 2 * numberOfMasses_; ++i)
	{
		out[i] = 0.0;
	}

	for (typeInt k = 1; k < numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}

		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k - 1] = p[0] / pSys_[0];
		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k] = - 2 * p[0] / pSys_[0];
		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k + 1] = p[0] / pSys_[0];

		for(typeInt i = numberOfMasses_ + k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	for(typeInt i = 0; i < 2 * numberOfMasses_ - 2; ++i)
	{
		out[(numberOfMasses_ - 1) * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[numberOfMasses_ * (2 * numberOfMasses_)  - 2] = p[0] / pSys_[0];
	out[numberOfMasses_ * (2 * numberOfMasses_)  - 1] = - 2 * p[0] / pSys_[0];

	for(typeInt i = 1; i < numberOfMasses_; ++i)
	{
		out[numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	for(typeInt i = numberOfMasses_ + 2; i < 2*numberOfMasses_; ++i)
	{
		out[numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[numberOfMasses_ * (2 * numberOfMasses_)] = 1.0;
	out[numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_] = - 2 * p[1] / pSys_[0];
	out[numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_ + 1] = p[1] / pSys_[0];

	for (typeInt k = numberOfMasses_ + 1; k < 2* numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
		out[k * (2 * numberOfMasses_) + k - numberOfMasses_] = 1.0;

		out[k * (2 * numberOfMasses_) + k - 1] = p[1] / pSys_[0];
		out[k * (2 * numberOfMasses_) + k] = - 2 * p[1] / pSys_[0];
		out[k * (2 * numberOfMasses_) + k + 1] = p[1] / pSys_[0];

		for(typeInt i = k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + numberOfMasses_ - 1] = 1.0;
	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 2] = p[1] / pSys_[0];
	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 1] = - 2 * p[1] / pSys_[0];
}

void MassSpringDamperProblemDescription::dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
	for(typeInt i = 0; i < numberOfMasses_; ++i)
	{
		out[i] = 0.0;
	}

	out[numberOfMasses_] = -2.0 / pSys_[0] * x[0] + 1.0 / pSys_[0] * x[1];

	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[numberOfMasses_ + k] = 1.0 / pSys_[0] * x[k - 1] - 2 / pSys_[0] * x[k] + 1.0 / pSys_[0] * x[k + 1];
	}

	out[2 * numberOfMasses_ - 1] = 1.0 / pSys_[0] * x[numberOfMasses_ - 2] - 2 / pSys_[0] * x[numberOfMasses_ - 1];

	for(typeInt i = 0; i < numberOfMasses_; ++i)
	{
		out[i + 2*numberOfMasses_] = 0.0;
	}

	out[3*numberOfMasses_] = - 2 / pSys_[0] * x[numberOfMasses_] + 1.0 / pSys_[0] * x[numberOfMasses_ + 1];

	for (typeInt k = 1; k <= numberOfMasses_ - 2; k++)
	{
		out[3*numberOfMasses_ + k] = 1.0 / pSys_[0] * x[numberOfMasses_ + k - 1] - 2 / pSys_[0] * x[numberOfMasses_ + k] + 1.0 / pSys_[0] * x[numberOfMasses_ + k + 1];
	}

	out[4 * numberOfMasses_ - 1] = 1.0 / pSys_[0] * x[2 * numberOfMasses_ - 2] - 2 / pSys_[0] * x[2 * numberOfMasses_ - 1];
}

void MassSpringDamperProblemDescription::dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p)
{
	for(typeInt i = 0; i < numberOfMasses_; ++i)
	{
		out[i] = 0.0;
	}

	out[numberOfMasses_] = -2 / pSys_[0];
	out[numberOfMasses_ + 1] = 1.0 / pSys_[0]; 

	for(typeInt i = numberOfMasses_ + 2; i < 2 * numberOfMasses_; ++i)
	{
		out[i] = 0.0;
	}

	for (typeInt k = 1; k < numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}

		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k - 1] = 1.0 / pSys_[0];
		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k] = - 2 / pSys_[0];
		out[k * (2 * numberOfMasses_) + numberOfMasses_ + k + 1] = 1.0 / pSys_[0];

		for(typeInt i = numberOfMasses_ + k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	for(typeInt i = 0; i < 2 * numberOfMasses_ - 2; ++i)
	{
		out[(numberOfMasses_ - 1) * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[numberOfMasses_ * (2 * numberOfMasses_)  - 2] = 1.0 / pSys_[0];
	out[numberOfMasses_ * (2 * numberOfMasses_)  - 1] = - 2 / pSys_[0];

	for(typeInt i = 1; i < numberOfMasses_; ++i)
	{
		out[numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	for(typeInt i = numberOfMasses_ + 2; i < 2*numberOfMasses_; ++i)
	{
		out[numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[numberOfMasses_ * (2 * numberOfMasses_)] = 1.0;
	out[numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_] = 0.0;
	out[numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_ + 1] = 0.0;

	for (typeInt k = numberOfMasses_ + 1; k < 2* numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
		out[k * (2 * numberOfMasses_) + k - numberOfMasses_] = 1.0;

		out[k * (2 * numberOfMasses_) + k - 1] = 0.0;
		out[k * (2 * numberOfMasses_) + k] = 0.0;
		out[k * (2 * numberOfMasses_) + k + 1] = 0.0;

		for(typeInt i = k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + numberOfMasses_ - 1] = 0.0;
	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 2] = 0.0;
	out[(2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 1] = 0.0;






	for(typeInt i = 0; i < numberOfMasses_; ++i)
	{
		out[4 * numberOfMasses_* numberOfMasses_ + i] = 0.0;
	}

	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_] = 0.0;
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ + 1] = 0.0; 

	for(typeInt i = numberOfMasses_ + 2; i < 2 * numberOfMasses_; ++i)
	{
		out[4 * numberOfMasses_* numberOfMasses_ + i] = 0.0;
	}

	for (typeInt k = 1; k < numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[4 * numberOfMasses_* numberOfMasses_ +  i + k * (2 * numberOfMasses_)] = 0.0;
		}

		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + numberOfMasses_ + k - 1] = 0.0;
		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + numberOfMasses_ + k] = 0.0;
		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + numberOfMasses_ + k + 1] = 0.0;

		for(typeInt i = numberOfMasses_ + k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[4 * numberOfMasses_* numberOfMasses_ + i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	for(typeInt i = 0; i < 2 * numberOfMasses_ - 2; ++i)
	{
		out[4 * numberOfMasses_* numberOfMasses_ + (numberOfMasses_ - 1) * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_)  - 2] = 0.0;
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_)  - 1] = 0.0;

	for(typeInt i = 1; i < numberOfMasses_; ++i)
	{
		out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	for(typeInt i = numberOfMasses_ + 2; i < 2*numberOfMasses_; ++i)
	{
		out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_) + i] = 0.0;
	}
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_)] = 0.0;
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_] = - 2 / pSys_[0];
	out[4 * numberOfMasses_* numberOfMasses_ + numberOfMasses_ * (2 * numberOfMasses_) + numberOfMasses_ + 1] = 1.0 / pSys_[0];

	for (typeInt k = numberOfMasses_ + 1; k < 2* numberOfMasses_ - 1; k++)
	{
		for(typeInt i = 0; i < numberOfMasses_ + k - 1; ++i)
		{
			out[4 * numberOfMasses_* numberOfMasses_ + i + k * (2 * numberOfMasses_)] = 0.0;
		}
		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + k - numberOfMasses_] = 0.0;

		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + k - 1] = 1.0 / pSys_[0];
		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + k] = - 2 / pSys_[0];
		out[4 * numberOfMasses_* numberOfMasses_ + k * (2 * numberOfMasses_) + k + 1] = 1.0 / pSys_[0];

		for(typeInt i = k + 2; i < 2 * numberOfMasses_; ++i)
		{
			out[4 * numberOfMasses_* numberOfMasses_ + i + k * (2 * numberOfMasses_)] = 0.0;
		}
	}

	out[4 * numberOfMasses_* numberOfMasses_ + (2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + numberOfMasses_ - 1] = 0.0;
	out[4 * numberOfMasses_* numberOfMasses_ + (2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 2] = 1.0 / pSys_[0];
	out[4 * numberOfMasses_* numberOfMasses_ + (2* numberOfMasses_ - 1) * (2 * numberOfMasses_) + 2 * numberOfMasses_ - 1] = - 2 / pSys_[0];
}