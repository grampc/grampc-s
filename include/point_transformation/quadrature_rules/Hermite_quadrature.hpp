#ifndef HERMITE_QUADRATURE_HPP
#define HERMITE_QUADRATURE_HPP

#include "polynomial/Hermite_polynomial_generator.hpp"
#include "quadrature_rule.hpp"


namespace grampc
{
    // Gauss–Hermite quadrature for a univariate Gaussian distribution, the measures of the Hermite polynomials is 1/sqrt(2*pi) * exp(-x^2 / 2)
    class HermiteQuadrature : public QuadratureRule
    {
    public:
        // Constructor for a quadrature of specified order
        HermiteQuadrature(typeInt quadratureOrder);

        // Roots of the corresponding orthogonal polynomials
         virtual const Vector& polynomialRoots() const override;

        // Get points for a normalized distribution with zero mean and variance = 1
        virtual const Vector& pointsNormalized() const override;
        
        // Get weights of the quadrature points 
        virtual const Vector& weights() const override;

        // Number of quadrature points
        virtual typeInt numberOfPoints() const override;

    private:
        // Points and weights for a univariate distribution 
        Vector roots_;
        Vector weights_;

        // Number of quadrature points
        typeInt quadratureOrder_;   
    };   
}

#endif // HERMITE_QUADRATURE_HPP