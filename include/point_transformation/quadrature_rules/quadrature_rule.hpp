#ifndef QUADRATURE_RULE_HPP
#define QUADRATURE_RULE_HPP

#include "util/grampc_s_constants.hpp"


namespace grampc
{
    // Univariate quadrature rule
    class QuadratureRule
    {
    public:
        virtual ~QuadratureRule() {};

        // Roots of the corresponding orthogonal polynomials
        virtual const Vector& polynomialRoots() const = 0;

        // Get points for a normalized distribution with zero mean and variance = 1
        virtual const Vector& pointsNormalized() const = 0;
        
        // Get weights of the quadrature points 
        virtual const Vector& weights() const = 0;

        // Number of quadrature points
        virtual typeInt numberOfPoints() const = 0;
    };   

    // Alias
    typedef std::shared_ptr<QuadratureRule> QuadratureRulePtr;
    typedef std::shared_ptr<const QuadratureRule> QuadratureRuleConstPtr;
}

#endif // QUADRATURE_RULE_HPP