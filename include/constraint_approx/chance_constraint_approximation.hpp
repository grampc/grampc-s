#ifndef CHANCE_CONSTRAINT_APPROXIMATION
#define CHANCE_CONSTRAINT_APPROXIMATION

#include <vector>
#include <memory>
#include "problem_description/grampc_interface.hpp"
#include "util/grampc_s_constants.hpp"


namespace grampc
{
    // Approximation of a chance constraint of the form z * sqrt(Var{h}) + E{h}
    class ChanceConstraintApproximation
    {
    public:
        virtual ~ChanceConstraintApproximation() {};

        // Set vector of chance constraint probabilities
        virtual void setConstraintProbability(const Vector& probabilities) = 0;

        // Get vector of chance constraint probabilities
        virtual const Vector& constraintProbability() const = 0;

        // Get coefficient z for the constraint tightening of the form z * sqrt(Var{h}) + E{h}
        virtual const Vector& tighteningCoefficient() const = 0;
    };

    // Alias
    typedef std::shared_ptr<ChanceConstraintApproximation> ChanceConstraintApproximationPtr;
    typedef std::shared_ptr<const ChanceConstraintApproximation> ChanceConstraintApproximationConstPtr;
}

#endif // CHANCE_CONSTRAINT_APPROXIMATION