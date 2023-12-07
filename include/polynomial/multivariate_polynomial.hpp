#ifndef MULTIVARIATE_POLYNOMIAL_HPP
#define MULTIVARIATE_POLYNOMIAL_HPP

#include "util/grampc_s_constants.hpp"
#include "polynomial.hpp"


namespace grampc
{
    // Multivariate polynomial defined by a product of univariate polynomials. 
    class MultivariatePolynomial
    {
    public:
        // Constructor for a multivariate polynomial
        MultivariatePolynomial(const std::vector<PolynomialConstPtr>& univariatePolynomials, const std::vector<typeRNum>& univariateSquaredNorm);

        // Evaluate the polynomial at a certain point
        typeRNum evaluate(const Vector& evaluationPoint) const;

        // Set all univariate polynomials
        void setPolynomials(std::vector<PolynomialConstPtr>& univariatePolynomials, std::vector<typeRNum>& univariateSquaredNorm);

        // Get polynomials
        const std::vector<PolynomialConstPtr>& polynomials() const;

        // Get one polynomial
        PolynomialConstPtr polynomials(typeInt index) const;

        // Get number of variables
        typeInt numVariables() const;

        // Get squared norm of the multivariate polynomial
        typeRNum squaredNorm() const;


    protected:
        std::vector<PolynomialConstPtr> polynomials_;
        typeInt numVariables_;
        typeRNum squaredNorm_;
    };

    // Alias
    typedef std::shared_ptr<MultivariatePolynomial> MultivariatePolynomialPtr;
    typedef std::shared_ptr<const MultivariatePolynomial> MultivariatePolynomialConstPtr;

    MultivariatePolynomialPtr MultivarPoly(const std::vector<PolynomialConstPtr>& univariatePolynomials, const std::vector<typeRNum>& univariateSquaredNorm);
}

#endif // MULTIVARIATE_POLYNOMIAL_HPP
