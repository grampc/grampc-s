#include "python_bindings.hpp"
#include "gaussian_process/stationary_kernel.hpp"
#include "gaussian_process/gaussian_process.hpp"
#include "gaussian_process/squared_exponential_kernel.hpp"
#include "gaussian_process/periodic_kernel.hpp"
#include "gaussian_process/locally_periodic_kernel.hpp"
#include "gaussian_process/Matern32_kernel.hpp"
#include "gaussian_process/Matern52_kernel.hpp"
#include "gaussian_process/kernel_sum.hpp"
#include "gaussian_process/kernel_product.hpp"

using namespace grampc;

template<typename Class>
py::class_<Class, StationaryKernel, std::shared_ptr<Class>> create_kernel(py::module_ &m, const std::string class_name)
{
    py::class_<Class, StationaryKernel, std::shared_ptr<Class>> cls(m, class_name.c_str());
    cls.def_property_readonly("dim", &Class::inputDimension)
        .def("evaluate", &Class::evaluate)
        .def("gradient", &Class::gradient);
    
    return cls;
}

void init_kernel(py::module_ &m)
{
    py::class_<StationaryKernel, StationaryKernelPtr>(m, "StationaryKernel");

    auto se = create_kernel<SquaredExponentialKernel>(m, "SquaredExponentialKernel");
    se.def(py::init<typeInt, typeRNum, const std::vector<typeRNum>&>());

    auto mat32 = create_kernel<Matern32Kernel>(m, "Matern32Kernel");
    mat32.def(py::init<typeInt, typeRNum, const std::vector<typeRNum>&>());
    
    auto mat52 = create_kernel<Matern52Kernel>(m, "Matern52Kernel");
    mat52.def(py::init<typeInt, typeRNum, const std::vector<typeRNum>&>());

    auto periodic = create_kernel<PeriodicKernel>(m, "PeriodicKernel");
    periodic.def(py::init<typeInt, typeRNum, const std::vector<typeRNum>&, const std::vector<typeRNum>&>());

    auto periodic_local = create_kernel<LocallyPeriodicKernel>(m, "LocallyPeriodicKernel");
    periodic_local.def(py::init<typeInt, typeRNum, const std::vector<typeRNum>&, const std::vector<typeRNum>&>());

    auto sum = create_kernel<KernelSum>(m, "KernelSum");
    sum.def(py::init<const std::vector<StationaryKernelPtr>&>());
    
    auto prod = create_kernel<KernelProduct>(m, "KernelProduct");
    prod.def(py::init<const StationaryKernelPtr&, const StationaryKernelPtr&>());
}

void init_gaussian_process(py::module_ &m)
{
    py::module_ m_kernel = m.def_submodule("kernel");
    init_kernel(m_kernel);

    py::class_<GaussianProcessData>(m, "GaussianProcessData")
        .def_readwrite("dim", &GaussianProcessData::inputDimension)
        .def_readwrite("size", &GaussianProcessData::numberOfDataPoints)
        .def_readwrite("output_variance", &GaussianProcessData::outputNoiseVariance)
        .def_readwrite("X", &GaussianProcessData::inputData)
        .def_readwrite("Y", &GaussianProcessData::outputData);

    py::class_<GaussianProcess, GaussianProcessPtr>(m, "GaussianProcess")
        .def(py::init<GaussianProcessData&, StationaryKernelPtr, const std::vector<bool>&, const std::vector<bool>&>())
        .def(py::init<const std::string&, StationaryKernelPtr, const std::vector<bool>&, const std::vector<bool>&>())
        .def("mean", &GaussianProcess::mean)
        .def("variance", &GaussianProcess::variance)
        .def("dmean_dx_vec", &GaussianProcess::dmean_dx_vec)
        .def("dmean_du_vec", &GaussianProcess::dmean_du_vec)
        .def("dvar_dx_vec", &GaussianProcess::dvar_dx_vec)
        .def("dvar_du_vec", &GaussianProcess::dvar_du_vec);
}