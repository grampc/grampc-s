#ifndef POLYNOMIAL_GENERATOR_HPP
#define POLYNOMIAL_GENERATOR_HPP

#include "polynomial.hpp"


namespace grampc
{
    // Generator for a family of orthogonal univariate polynomials
    class OrthogonalPolynomialGenerator
    {
    public:
        virtual ~OrthogonalPolynomialGenerator() {};

        // Get a polynomial
        virtual PolynomialConstPtr getPolynomial(typeInt order) const = 0;

        // Get the squared norm of the polynomial 
        virtual typeRNum getSquaredNorm(typeInt order) const = 0;

        // Get maximum order of polynomial that can be created
        virtual typeInt getMaximumOrder() const = 0;
    };   

    // Alias
    typedef std::shared_ptr<OrthogonalPolynomialGenerator> PolynomialGeneratorPtr;
}

#endif // POLYNOMIAL_GENERATOR_HPP