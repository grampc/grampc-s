#ifndef GRAMPC_S_UTIL_HPP
#define GRAMPC_S_UTIL_HPP

#include <vector>
#include <iostream>
#include "grampc_s_constants.hpp"
#include "point_transformation/quadrature_rules/Hermite_quadrature.hpp"
#include "point_transformation/quadrature_rules/Legendre_quadrature.hpp"
#include "polynomial/Hermite_polynomial_generator.hpp"
#include "polynomial/Legendre_polynomial_generator.hpp"

namespace grampc
{
    // Factorial
    long long factorial(long long n);

    // Return the quadrature rule that corresponds to the distribution
    QuadratureRuleConstPtr correspondingQuadratureRule(PolynomialFamily polyFam, typeInt quadratureOrder);

    // Return the quadrature rule that corresponds to the distribution
    PolynomialGeneratorPtr correspondingPolynomialGenerator(PolynomialFamily polyFam, typeInt maxOrder);

    void deriveCholesky(VectorRef out, MatrixConstRef chol, typeInt k, typeInt l);
  }

#endif // GRAMPC_S_UTIL_HPP
