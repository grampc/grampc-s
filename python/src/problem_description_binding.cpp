#include "python_bindings.hpp"
#include "problem_description/resampling_problem_description.hpp"
#include "problem_description/monte_carlo_problem_description.hpp"
#include "problem_description/resampling_GP_problem_description.hpp"
#include "problem_description/sigma_point_problem_description.hpp"
#include "problem_description/taylor_problem_description.hpp"

using namespace grampc;

void init_problem_descriptions(py::module_ &m)
{
    py::module_::import("pygrampc");
    
    py::class_<TaylorBaseProblemDescription, ProblemDescription, TaylorBaseProblemDescriptionPtr> taylorbase(m, "TaylorBaseProblemDescription");
    
    py::class_<ResamplingProblemDescription, ProblemDescription, ResamplingProblemDescriptionPtr> resampling(m, "ResamplingProblemDescription");
    resampling.def(py::init<ProblemDescriptionPtr, ChanceConstraintApproximationConstPtr, PointTransformationPtr>())
        .def(py::init<ProblemDescriptionPtr, ChanceConstraintApproximationConstPtr, PointTransformationPtr, MatrixConstRef>())
        .def(py::init<ProblemDescriptionPtr, PointTransformationPtr>())
        .def(py::init<ProblemDescriptionPtr, PointTransformationPtr, MatrixConstRef>());

    problem_methods<ResamplingProblemDescription>(resampling);

    py::class_<ResamplingGPProblemDescription, ProblemDescription, ResamplingGPProblemDescriptionPtr> resampling_gp(m, "ResamplingGPProblemDescription");
    resampling_gp.def(py::init<ProblemDescriptionPtr, ChanceConstraintApproximationConstPtr, PointTransformationPtr, const std::vector<GaussianProcessPtr>&, const std::vector<typeInt>&>())
        .def(py::init<ProblemDescriptionPtr, PointTransformationPtr,const std::vector<GaussianProcessPtr>&, const std::vector<typeInt>&>());

    problem_methods<ResamplingGPProblemDescription>(resampling_gp);
        
    py::class_<MonteCarloProblemDescription, ProblemDescription, MonteCarloProblemDescriptionPtr> monte_carlo(m, "MonteCarloProblemDescription");
    monte_carlo.def(py::init<ProblemDescriptionPtr, PointTransformationPtr>());

    problem_methods<MonteCarloProblemDescription>(monte_carlo);

    py::class_<SigmaPointProblemDescription, ProblemDescription, SigmaPointProblemDescriptionPtr> sigma_point(m, "SigmaPointProblemDescription");
    sigma_point.def(py::init<ProblemDescriptionPtr, ChanceConstraintApproximationConstPtr, PointTransformationPtr>())
        .def(py::init<ProblemDescriptionPtr, PointTransformationPtr>());

    problem_methods<SigmaPointProblemDescription>(sigma_point);

    py::class_<TaylorProblemDescription, ProblemDescription, TaylorProblemDescriptionPtr> taylor(m, "TaylorProblemDescription");
    taylor.def(py::init<TaylorBaseProblemDescriptionPtr, ChanceConstraintApproximationConstPtr>())
        .def(py::init<TaylorBaseProblemDescriptionPtr, ChanceConstraintApproximationConstPtr, MatrixConstRef>())
        .def(py::init<TaylorBaseProblemDescriptionPtr>())
        .def(py::init<TaylorBaseProblemDescriptionPtr, MatrixConstRef>());

    problem_methods<TaylorProblemDescription>(taylor);
}