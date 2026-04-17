#include "python_bindings.hpp"
#include "constraint_approx/chance_constraint_approximation.hpp"
#include "constraint_approx/Gaussian_constraint_approximation.hpp"
#include "constraint_approx/chebyshev_constraint_approximation.hpp"
#include "constraint_approx/symmetric_constraint_approximation.hpp"

using namespace grampc;

template<typename Class>
void create_constraint(py::module_ &m, const std::string class_name)
{
    py::class_<Class, ChanceConstraintApproximation, std::shared_ptr<Class>> cls(m, class_name.c_str());
    cls.def(py::init<const Vector&>())
        .def_property("probability", &Class::setConstraintProbability, &Class::constraintProbability)
        .def_property_readonly("z", &Class::tighteningCoefficient);
}

void init_constraints(py::module_ &m)
{
    py::class_<ChanceConstraintApproximation, ChanceConstraintApproximationPtr>(m, "ChanceConstraintApproximation");

    create_constraint<GaussianConstraintApproximation>(m, "GaussianConstraintApproximation");
    create_constraint<ChebyshevConstraintApproximation>(m, "ChebyshevConstraintApproximation");
    create_constraint<SymmetricConstraintApproximation>(m, "SymmetricConstraintApproximation");
}