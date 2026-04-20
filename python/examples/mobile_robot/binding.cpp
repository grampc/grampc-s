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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "mobile_robot_problem_description.hpp"
#include "mobile_robot_with_parameters_problem_description.hpp"

namespace py = pybind11;

PYBIND11_MODULE(grampc_s_mobile_robot, m)
{
    py::module_::import("pygrampc");

    py::class_<MobileRobotProblemDescription, grampc::ProblemDescription, std::shared_ptr<MobileRobotProblemDescription>> mobile_robot_problem(m, "MobileRobotProblem");
    mobile_robot_problem.def(py::init<const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::vector<typeRNum>&>());

    py::class_<MobileRobotWithParametersProblemDescription, grampc::ProblemDescription, std::shared_ptr<MobileRobotWithParametersProblemDescription>> mobile_robot_with_parameters_problem(m, "MobileRobotWithParametersProblem");
    mobile_robot_with_parameters_problem.def(py::init<const std::vector<typeRNum>&, const std::vector<typeRNum>&, const std::vector<typeRNum>&>());
}