#include "python_bindings.hpp"
#include "distribution/distribution.hpp"
#include "distribution/Gaussian_distribution.hpp"
#include "distribution/univariate_chi_squared_distribution.hpp"
#include "distribution/univariate_exponential_distribution.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"
#include "distribution/univariate_piecewise_constant_distribution.hpp"
#include "distribution/univariate_uniform_distribution.hpp"
#include "distribution/univariate_gamma_distribution.hpp"
#include "distribution/univariate_lognormal_distribution.hpp"
#include "distribution/univariate_student-t_distribution.hpp"
#include "distribution/univariate_weibull_distribution.hpp"
#include "distribution/univariate_extreme_value_distribution.hpp"
#include "distribution/univariate_f_distribution.hpp"
#include "distribution/univariate_beta_distribution.hpp"

using namespace grampc;

template<typename Class>
py::class_<Class, Distribution, std::shared_ptr<Class>> create_distribution(py::module_ &m, const std::string class_name)
{
    py::class_<Class, Distribution, std::shared_ptr<Class>> cls(m, class_name.c_str());
    cls.def_property_readonly("dimension", &Class::dimension)
        .def_property_readonly("mean", &Class::mean)
        .def_property_readonly("covariance", &Class::covariance)
        .def_property_readonly("cov_cholesky", &Class::covCholesky)
        .def_property_readonly("polynomial_families", &Class::polynomialFamily)
        .def("sample", &Class::sample)
        .def("set_mean_and_covariance", &Class::setMeanAndCovariance);
    return cls;
}

void init_distributions(py::module_ &m)
{
    py::class_<Distribution, DistributionPtr>(m, "Distribution")
        .def(py::init<typeInt>())
        .def(py::init<VectorConstRef, MatrixConstRef, const std::vector<PolynomialFamily>&>())
        .def(py::init<VectorConstRef, MatrixConstRef, MatrixConstRef, const std::vector<PolynomialFamily>&>())
        .def_property_readonly("dimension", &Distribution::dimension)
        .def_property_readonly("mean", &Distribution::mean)
        .def_property_readonly("covariance", &Distribution::covariance)
        .def_property_readonly("cov_cholesky", &Distribution::covCholesky)
        .def_property_readonly("polynomial_families", &Distribution::polynomialFamily)
        .def("sample", &Distribution::sample)
        .def("set_mean_and_covariance", &Distribution::setMeanAndCovariance);

    auto multivariate = create_distribution<MultivariateDistribution>(m, "MultivariateDistribution");
    multivariate.def(py::init<const std::vector<DistributionPtr>>())
        .def("replace_distribution", &MultivariateDistribution::replaceDistribution)
        .def("get_distribution", &MultivariateDistribution::getDistribution);

    auto gauss = create_distribution<GaussianDistribution>(m, "GaussianDistribution");
    gauss.def(py::init<VectorConstRef, MatrixConstRef>())
        .def(py::init<typeRNum, typeRNum>());

    auto uni_beta = create_distribution<UnivariateBetaDistribution>(m, "UnivariateBetaDistribution");
    uni_beta.def(py::init<typeRNum, typeRNum>());

    auto uni_chi = create_distribution<UnivariateChiSquaredDistribution>(m, "UnivariateChiSquaredDistribution");
    uni_chi.def(py::init<typeRNum>());

    auto uni_exp = create_distribution<UnivariateExponentialDistribution>(m, "UnivariateExponentialDistribution");
    uni_exp.def(py::init<typeRNum>());

    auto uni_extr = create_distribution<UnivariateExtremeValueDistribution>(m, "UnivariateExtremeValueDistribution");
    uni_extr.def(py::init<typeRNum, typeRNum>());

    auto uni_f = create_distribution<UnivariateFDistribution>(m, "UnivariateFDistribution");
    uni_f.def(py::init<typeRNum, typeRNum>());

    auto uni_gamma = create_distribution<UnivariateGammaDistribution>(m, "UnivariateGammaDistribution");
    uni_gamma.def(py::init<typeRNum, typeRNum>());

    auto uni_log_normal = create_distribution<UnivariateLognormalDistribution>(m, "UnivariateLogNormalDistribution");
    uni_log_normal.def(py::init<typeRNum, typeRNum>());

    auto uni_piece = create_distribution<UnivariatePiecewiseConstantDistribution>(m, "UnivariatePiecewiseConstantDistribution");
    uni_piece.def(py::init<const std::vector<typeRNum>&, const std::vector<typeRNum>&>());

    auto uni_student = create_distribution<UnivariateStudentTDistribution>(m, "UnivariateStudentTDistribution");
    uni_student.def(py::init<typeRNum, typeRNum, typeRNum>(), py::arg("nu"), py::arg("mu") = 0.0, py::arg("sigma") = 1.0);

    auto uni_uniform = create_distribution<UnivariateUniformDistribution>(m, "UnivariateUniformDistribution");
    uni_uniform.def(py::init<typeRNum, typeRNum>());

    auto uni_weibull = create_distribution<UnivariateWeibullDistribution>(m, "UnivariateWeibullDistribution");
    uni_weibull.def(py::init<typeRNum, typeRNum>());
}
