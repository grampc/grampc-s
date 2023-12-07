#ifndef LEGENDRE_POLYNOMIAL_GENERATOR_HPP
#define LEGENDRE_POLYNOMIAL_GENERATOR_HPP

#include <algorithm>
#include "orthogonal_polynomial_generator.hpp"


namespace grampc
{
    // Generator for univariate Legendre polynomials with measure 1/2
    class LegendrePolynomialGenerator : public OrthogonalPolynomialGenerator
    {
    public:
        // Constructor specifying the maximum polynomial order
        LegendrePolynomialGenerator(typeInt maxOrder);

        // Get a polynomial
        virtual PolynomialConstPtr getPolynomial(typeInt order) const override;

        // Get the squared norm of the polynomial 
        virtual typeRNum getSquaredNorm(typeInt order) const override;

        // Get maximum order of polynomial that can be created
        virtual typeInt getMaximumOrder() const override;

    private:
        std::vector<PolynomialPtr> legendrePolynomials;
        Vector legendreSquaredNorm;
    };   
}


#endif //LEGENDRE_POLYNOMIAL_GENERATOR_HPP