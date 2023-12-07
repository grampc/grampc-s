#ifndef LEGENDRE_QUADRATURE_POINT_GENERATOR_HPP
#define LEGENDRE_QUADRATURE_POINT_GENERATOR_HPP

#include "polynomial/Legendre_polynomial_generator.hpp"
#include "quadrature_rule.hpp"


namespace grampc
{
    // Gauss–Legendre quadrature for a univariate uniform distribution, the measures of the Legendre polynomials is 1/2
    class LegendreQuadrature : public QuadratureRule
    {
    public:
        // Constructor for a quadrature of specified order
        LegendreQuadrature(typeInt quadratureOrder);

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
        Vector pointsNormalized_;
        Vector roots_;
        Vector weights_;

        // Number of quadrature points
        typeInt quadratureOrder_;       
    };   
}

#endif // LEGENDRE_QUADRATURE_POINT_GENERATOR_HPP