# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
    grampc/src/discrete.c
    grampc/src/finite_diff.c
	grampc/src/grampc_erk.c
	grampc/src/grampc_alloc.c
    grampc/src/grampc_configreader.c
	grampc/src/grampc_fixedsize.c
	grampc/src/grampc_init.c
	grampc/src/grampc_mess.c
	grampc/src/grampc_run.c
	grampc/src/grampc_setopt.c
	grampc/src/grampc_setparam.c
	grampc/src/grampc_util.c
	grampc/src/rodas.c
	grampc/src/ruku45.c
	grampc/src/simpson.c
	grampc/src/trapezoidal.c)