#include "python_bindings.hpp"
#include "simulator/simulator.hpp"
#include "util/grampc_s_constants.hpp"

using namespace grampc;

PYBIND11_MODULE(grampc_s, m)
{
    py::module_ m_distributions = m.def_submodule("distributions");
    py::module_ m_constraints = m.def_submodule("constraints");
    py::module_ m_transforms = m.def_submodule("transformations");
    py::module_ m_polynomial = m.def_submodule("polynomials");
    py::module_ m_problem_descriptions = m.def_submodule("problem_descriptions");
    py::module_ m_gaussian_process = m.def_submodule("gaussian_process");
    
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<VectorConstRef, VectorConstRef, typeInt, const SystemFct, const std::string&, MatrixConstRef, typeRNum, typeRNum, typeRNum, bool>())
        .def(py::init<VectorConstRef, typeInt, const SystemFct, const std::string&, MatrixConstRef, typeRNum, typeRNum, typeRNum, bool>())
        .def(py::init<VectorConstRef, VectorConstRef, typeInt, const SystemFct, const std::string&, typeRNum, typeRNum, typeRNum, bool>())
        .def(py::init<VectorConstRef, typeInt, const SystemFct, const std::string&, typeRNum, typeRNum, typeRNum, bool>())
        .def("simulate", &Simulator::simulate)
        .def_property_readonly("t", &Simulator::time);

    py::class_<RandomNumberGenerator>(m, "RandomNumberGenerator")
        .def(py::init<typeInt>());

    init_polynomial(m_polynomial);
    init_distributions(m_distributions);
    init_constraints(m_constraints);
    init_transformations(m_transforms);
    init_gaussian_process(m_gaussian_process);
    init_problem_descriptions(m_problem_descriptions);
}
