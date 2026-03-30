#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/native_enum.h>
#include "problem_description/problem_description.hpp"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void init_distributions(py::module_ &);
void init_constraints(py::module_ &);
void init_transformations(py::module_ &);
void init_problem_descriptions(py::module_ &);
void init_gaussian_process(py::module_ &);
void init_kernel(py::module_ &);
void init_polynomial(py::module &);

template <typename Class>
void problem_methods(py::class_<Class, grampc::ProblemDescription, std::shared_ptr<Class>> &cls)
{
    cls.def("compute_x0_and_p0", overload_cast_<DistributionPtr>()(&Class::compute_x0_and_p0));
    cls.def("compute_x0_and_p0", overload_cast_<DistributionPtr, DistributionPtr>()(&Class::compute_x0_and_p0));
    cls.def_property_readonly("x0", &Class::x0);
    cls.def_property_readonly("p0", &Class::p0);
}
