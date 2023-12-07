#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <vector>
#include <memory>
#include "util/grampc_s_constants.hpp"



namespace grampc
{
    // A univariate polynomial defined by its coefficients. The first coefficient defines the monomial of order 0, the last coefficient defines the monomial of the highest order
    class Polynomial
    {
    public:
        // Constructor for an empty polynomial object
        Polynomial();

        // Constructor for an polynomial definde by its coefficients
        Polynomial(const Vector& coefficients);
        
        // Evaluate the polynomial at one or multiple points
        typeRNum evaluate(ctypeRNum evaluationPoint) const;

        // Evaluate the gradient of the polynomial at one or multiple points
        typeRNum gradient(ctypeRNum evaluationPoint) const;

        // Evaluate the hessian of the polynomial at one or multiple points
        typeRNum hessian(ctypeRNum evaluationPoint) const;

        // Add a polynomial
        void addPolynomial(const Polynomial& polynomial);

        // Subtract a polynomial
        void subtractPolynomial(const Polynomial& polynomial);

        // Multiply a polynomial
        void multiplyPolynomial(const Polynomial& polynomial);

        // Multiply with a scalar
        void multiplyScalar(const typeRNum factor);

        // Set coefficients
        void setCoefficients(const Vector& coefficients);

        // Get coefficients
        const Vector& getCoefficients() const;

        // Get coefficient
        typeRNum getCoefficient(int index) const;

        // Get number of coefficients
        typeInt getNumCoefficients() const;


    protected:
        typeInt numCoefficients_;
        Vector coefficients_;
    };

    // Alias
    typedef std::shared_ptr<Polynomial> PolynomialPtr;
    typedef std::shared_ptr<const Polynomial> PolynomialConstPtr;

    PolynomialPtr Poly();
    PolynomialPtr Poly(const Vector& coefficients);
}

#endif // POLYNOMIAL_HPP
