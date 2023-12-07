#ifndef HERMITE_POLYNOMIAL_GENERATOR_HPP
#define HERMITE_POLYNOMIAL_GENERATOR_HPP

#include <algorithm>
#include "orthogonal_polynomial_generator.hpp"
#include "util/grampc_s_util.hpp"


namespace grampc
{
    // Generator for univariate Hermite polynomials with measure 1/sqrt(2*pi) * exp(-x^2 / 2) 
    class HermitePolynomialGenerator : public OrthogonalPolynomialGenerator
    {
    public:
        // Constructor specifying the maximum polynomial order
        HermitePolynomialGenerator(typeInt maxOrder);

        // Get a polynomial
        virtual PolynomialConstPtr getPolynomial(typeInt order) const override;

        // Get the squared norm of the polynomial 
        virtual typeRNum getSquaredNorm(typeInt order) const override;

        // Get maximum order of polynomial that can be created
        virtual typeInt getMaximumOrder() const override;

    private:
        std::vector<PolynomialPtr> hermitePolynomials;
        Vector hermiteSquaredNorm;
    };   
}


#endif //HERMITE_POLYNOMIAL_GENERATOR_HPP