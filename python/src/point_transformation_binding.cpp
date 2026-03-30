#include "python_bindings.hpp"
#include "point_transformation/point_transformation.hpp"
#include "point_transformation/unscented_transformation.hpp"
#include "point_transformation/monte_carlo.hpp"
#include "point_transformation/PCE_transformation.hpp"
#include "point_transformation/stirling_interpolation_first_order.hpp"
#include "point_transformation/stirling_interpolation_second_order.hpp"
#include "point_transformation/composed_quadrature.hpp"

using namespace grampc;

template<typename Class>
void transform_methods(py::class_<Class, PointTransformation, std::shared_ptr<Class>> &cls)
{
    cls.def("mean", &Class::mean);
    cls.def("covariance", &Class::covariance);
    cls.def("variance", &Class::variance);
    cls.def("mean1D", &Class::mean1D);
    cls.def("dmean_dpoints_vec", &Class::dmean_dpoints_vec);
    cls.def("dcov_dpointsX_vec", &Class::dcov_dpointsX_vec);
    cls.def("dcov_dpointsY_vec", &Class::dcov_dpointsY_vec);
    cls.def("dpoints_dmean_vec", &Class::dpoints_dmean_vec);
    cls.def("dpoints_dcov_vec", &Class::dpoints_dcov_vec);
    cls.def("dmean1D_dpoints", &Class::dmean1D_dpoints);
    cls.def("dvar_points", &Class::dvar_dpoints);
    cls.def("number_of_points", &Class::numberOfPoints);
}

void init_transformations(py::module_ &m)
{
    py::class_<PointTransformation, PointTransformationPtr>(m, "PointTransformation");

    py::class_<UnscentedTransformation, PointTransformation, std::shared_ptr<UnscentedTransformation>> unscented_transform(m, "UnscentedTransformation");
    unscented_transform.def(py::init<typeInt, typeInt, typeRNum, typeRNum, typeRNum>())
        .def(py::init<typeInt, typeInt, typeRNum, typeRNum, typeRNum, const std::vector<bool>&>())
        .def("points", overload_cast_<VectorConstRef, MatrixConstRef>()(&UnscentedTransformation::points));
    
    transform_methods<UnscentedTransformation>(unscented_transform);
    
    py::class_<MonteCarloTransformation, PointTransformation, std::shared_ptr<MonteCarloTransformation>> monte_carlo_transform(m, "MonteCarloTransformation");
    monte_carlo_transform.def(py::init<typeInt, typeInt, typeInt, const RandomNumberGenerator&>())
        .def("points", overload_cast_<>()(&MonteCarloTransformation::points))
        .def("points", overload_cast_<DistributionConstPtr>()(&MonteCarloTransformation::points))
        .def("set_point", &MonteCarloTransformation::setPoint);

    transform_methods<MonteCarloTransformation>(monte_carlo_transform);
        
    py::class_<PCE_Transformation, PointTransformation, std::shared_ptr<PCE_Transformation>> pce(m, "PCE_Transformation");
    pce.def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, typeInt, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>&>())
        .def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, typeInt, const std::vector<typeInt>&>())
        .def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, typeInt, typeInt>())
        .def("points", overload_cast_<>()(&PCE_Transformation::points))
        .def("points", overload_cast_<DistributionConstPtr>()(&PCE_Transformation::points));

    transform_methods<PCE_Transformation>(pce);

    py::class_<StirlingInterpolationFirstOrder, PointTransformation, std::shared_ptr<StirlingInterpolationFirstOrder>> stirling_first(m, "StirlingInterpolationFirstOrder");
    stirling_first.def(py::init<typeInt, typeInt, typeRNum, const std::vector<bool>&>())
        .def(py::init<typeInt, typeInt, typeRNum>())
        .def("points", overload_cast_<>()(&StirlingInterpolationFirstOrder::points))
        .def("points", overload_cast_<DistributionConstPtr>()(&StirlingInterpolationFirstOrder::points))
        .def("points", overload_cast_<VectorConstRef, MatrixConstRef>()(&StirlingInterpolationFirstOrder::points));

    transform_methods<StirlingInterpolationFirstOrder>(stirling_first);
    
    py::class_<StirlingInterpolationSecondOrder, PointTransformation, std::shared_ptr<StirlingInterpolationSecondOrder>> stirling_second(m, "StirlingInterpolationSecondOrder");
    stirling_second.def(py::init<typeInt, typeInt, typeRNum, const std::vector<bool>&>())
        .def(py::init<typeInt, typeInt, typeRNum>())
        .def("points", overload_cast_<>()(&StirlingInterpolationSecondOrder::points))
        .def("points", overload_cast_<DistributionConstPtr>()(&StirlingInterpolationSecondOrder::points))
        .def("points", overload_cast_<VectorConstRef, MatrixConstRef>()(&StirlingInterpolationSecondOrder::points));

    transform_methods<StirlingInterpolationSecondOrder>(stirling_second);

    py::class_<ComposedQuadrature, PointTransformation, std::shared_ptr<ComposedQuadrature>> composed(m, "ComposedQuadrature");
    composed.def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, const Eigen::Ref<const Eigen::Vector<typeInt, Eigen::Dynamic>>&>())
        .def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, const std::vector<typeInt>&>())
        .def(py::init<typeInt, typeInt, const std::vector<PolynomialFamily>&, typeInt>())
        .def("points", overload_cast_<>()(&ComposedQuadrature::points))
        .def("points", overload_cast_<DistributionConstPtr>()(&ComposedQuadrature::points))
        .def("points", overload_cast_<VectorConstRef, MatrixConstRef>()(&ComposedQuadrature::points))
        .def("roots", &ComposedQuadrature::roots)
        .def("weights", &ComposedQuadrature::weights)
        .def("normalized_points", &ComposedQuadrature::normalizedPoints);

    transform_methods<ComposedQuadrature>(composed);
}