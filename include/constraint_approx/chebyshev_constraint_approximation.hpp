#ifndef CHEBYSHEV_CONSTRAINT_APPROXIMATION
#define CHEBYSHEV_CONSTRAINT_APPROXIMATION

#include "chance_constraint_approximation.hpp"


namespace grampc
{
    // Chance constraint approximation using the Chebyshev inequality
    class ChebyshevConstraintApproximation : public ChanceConstraintApproximation
    {
    public:
        // Constructor for a constraint vector
        ChebyshevConstraintApproximation(const Vector& probabilities);

        // Set vector of chance constraint probabilities
        virtual void setConstraintProbability(const Vector& probabilities) override;

        // Get vector of chance constraint probabilities
        virtual const Vector& constraintProbability() const override;

        // Get coefficient z for the constraint tightening of the form z * sqrt(Var{h}) + E{h}
        virtual const Vector& tighteningCoefficient() const override;

    private:
        Vector probabilities_;
        Vector coefficients_;
    };

    ChanceConstraintApproximationPtr Chebyshev(const Vector& probabilities);
}


#endif // CHEBYSHEV_CONSTRAINT_APPROXIMATION