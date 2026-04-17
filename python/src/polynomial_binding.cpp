#include "python_bindings.hpp"
#include "polynomial/polynomial.hpp"
#include "polynomial/multivariate_polynomial.hpp"
#include "polynomial/orthogonal_polynomial_generator.hpp"
#include "polynomial/Hermite_polynomial_generator.hpp"
#include "polynomial/Legendre_polynomial_generator.hpp"

using namespace grampc;

void init_polynomial(py::module_ &m)
{
    py::enum_<PolynomialFamily>(m, "PolynomialFamily")
        .value("NONE", PolynomialFamily::NONE)
        .value("HERMITE", PolynomialFamily::HERMITE)
        .value("LEGENDRE", PolynomialFamily::LEGENDRE)
        .export_values();

    py::class_<Polynomial, PolynomialPtr>(m, "Polynomial")
        .def(py::init<>())
        .def(py::init<const Vector&>())
        .def("evaluate", &Polynomial::evaluate)
        .def("gradient", &Polynomial::gradient)
        .def("hessian", &Polynomial::hessian)
        .def("add_polynomial", &Polynomial::addPolynomial)
        .def("subtract_polynomial", &Polynomial::subtractPolynomial)
        .def("multiply_polynomial", &Polynomial::multiplyPolynomial)
        .def("multiply_scalar", &Polynomial::multiplyScalar)
        .def_property("coefficients", &Polynomial::getCoefficients, &Polynomial::setCoefficients)
        .def_property_readonly("num_coefficients", &Polynomial::getNumCoefficients);

    py::class_<MultivariatePolynomial, MultivariatePolynomialPtr>(m, "MultivariatePolynomial")
        .def(py::init<const std::vector<PolynomialConstPtr>&, const std::vector<typeRNum>&>())
        .def("evaluate", &MultivariatePolynomial::evaluate)
        .def_property_readonly("polynomials", &MultivariatePolynomial::polynomials)
        .def_property_readonly("num_variables", &MultivariatePolynomial::numVariables)
        .def_property_readonly("squared_norm", &MultivariatePolynomial::squaredNorm);

    py::class_<OrthogonalPolynomialGenerator>(m, "OrthogonalPolynomialGenerator");

    py::class_<HermitePolynomialGenerator, OrthogonalPolynomialGenerator>(m, "HermitePolynomialGenerator")
        .def(py::init<typeInt>())
        .def("get_polynomial", &HermitePolynomialGenerator::getPolynomial)
        .def("get_squared_norm", &HermitePolynomialGenerator::getSquaredNorm)
        .def("get_maximum_order", &HermitePolynomialGenerator::getMaximumOrder);

    py::class_<LegendrePolynomialGenerator, OrthogonalPolynomialGenerator>(m, "LegendrePolynomialGenerator")
        .def(py::init<typeInt>())
        .def("get_polynomial", &LegendrePolynomialGenerator::getPolynomial)
        .def("get_squared_norm", &LegendrePolynomialGenerator::getSquaredNorm)
        .def("get_maximum_order", &LegendrePolynomialGenerator::getMaximumOrder);

}