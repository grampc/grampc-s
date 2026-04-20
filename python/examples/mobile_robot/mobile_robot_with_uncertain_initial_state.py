import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from grampc_s.constraints import GaussianConstraintApproximation
from grampc_s.distributions import GaussianDistribution
from grampc_s.transformations import UnscentedTransformation
from grampc_s.problem_descriptions import SigmaPointProblemDescription
from pygrampc import Grampc, GrampcResults
from grampc_s_mobile_robot import MobileRobotProblem

Nx = 3
Nc = 4

Parameters = {
    "umax": [2, 2],
    "umin": [-2,-2],
    "Thor": 1.0,
    "dt": 0.01
}

Options = {
    "Nhor": 30,
    "MaxGradIter": 3,
    "MaxMultIter": 1,
    "InequalityConstraints": "on"
}

if __name__ == "__main__":
    plotSteps = 10

    # problem-specific parameters
    pSys = []
    pCost = [1, 5, 0.1, 0.125, 0.0125]
    pCon = [-0.9, 0.9, -0.9, 0.9]

    # deterministic problem description
    problem = MobileRobotProblem(pSys, pCost, pCon)
    # approximation of chance constraints
    chance_constraint = GaussianConstraintApproximation([0.95] * Nc)
    # method for uncertainty propagation
    point_transformation = UnscentedTransformation(Nx, Nx, 1, 2, 0)
    # stochastic problem description
    stochastic_problem = SigmaPointProblemDescription(problem, chance_constraint, point_transformation)
    # create solver
    grampc = Grampc(stochastic_problem) # instead of Grampc(problem)

    # distribution of initial state
    mean = np.array([0, 0, 0])
    covariance = np.diag(np.array([1e-2, 1e-2, 1e-2]))
    initial_state = GaussianDistribution(mean, covariance)
    # compute initial state of approximation
    stochastic_problem.compute_x0_and_p0(initial_state)
    grampc.set_param({'x0': stochastic_problem.x0})

    # set parameters and options for GRAMPC
    grampc.set_param(Parameters)
    grampc.set_opt(Options)
    grampc.estim_penmin(True)

    # MPC main loop
    Tsim = 4.5
    vec = GrampcResults(grampc, Tsim)

    vec.mean = np.full((len(vec.t), Nx), np.nan)
    vec.std  = np.full((len(vec.t), Nx), np.nan)
    pred_mean = np.full((grampc.opt.Nhor, Nx), np.nan)
    pred_std = np.full((grampc.opt.Nhor, Nx), np.nan)
    zalpha = chance_constraint.z[0] # assuming all constraints have same threshold

    # prepare plots (interactive mode)
    plt.ion()
    fig, axs = plt.subplots(1, 2)
    axs[0].add_patch(plt.Circle((0, 0), 1, fill=False, linestyle='--'))
    axs[0].add_patch(plt.Rectangle((-0.9, -0.9), 1.8, 1.8, fill=False, linestyle=':'))
    ph_vec = axs[0].plot(vec.mean[:, 0], vec.mean[:, 1])[0]
    ph_pred = [0] * point_transformation.number_of_points()
    for j in range(point_transformation.number_of_points()):
        ph_pred[j] = axs[0].plot(grampc.rws.x[Nx*j, :], grampc.rws.x[Nx*j+1, :], '.')[0]
    axs[0].axis('equal')
    axs[0].set_xlabel('Position x')
    axs[0].set_ylabel('Position y')

    axs[1].plot([0, Tsim+grampc.param.Thor], [0.9, 0.9], 'k--')
    axs[1].plot([0, Tsim+grampc.param.Thor], [-0.9, -0.9], 'k--')
    ph_vec_mean = axs[1].plot(vec.t, vec.mean)
    axs[1].set_prop_cycle(None)
    ph_vec_std_up = axs[1].plot(vec.t, vec.mean + zalpha * vec.std, '--')
    axs[1].set_prop_cycle(None)
    ph_vec_std_down = axs[1].plot(vec.t, vec.mean - zalpha * vec.std, '--')
    axs[1].set_prop_cycle(None)
    ph_pred_mean = axs[1].plot(grampc.rws.t, pred_mean, '.')
    axs[1].set_prop_cycle(None)
    ph_pred_std_up = axs[1].plot(grampc.rws.t, pred_mean + zalpha * pred_std, ':')
    axs[1].set_prop_cycle(None)
    ph_pred_std_down = axs[1].plot(grampc.rws.t, pred_mean - zalpha * pred_std, ':')
    axs[1].set_prop_cycle(None)
    axs[1].set_xlim(0, Tsim+grampc.param.Thor)
    axs[1].set_ylim(-1.1, 1.1)
    axs[1].set_xlabel('Time t')
    axs[1].set_ylabel('States x')

    for i, t in enumerate(vec.t):
        # solve problem
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)
        vec.x[i, :] = grampc.param.x0 # store current state instead of predicted state

        # compute mean and standard deviation of simulated states
        points = np.reshape(vec.x[i, :], (Nx, point_transformation.number_of_points()), order='F')
        P = point_transformation.covariance(points, points)
        vec.mean[i, :] = point_transformation.mean(points)
        vec.std[i, :] = np.sqrt(np.diag(P))
        # compute mean and standard deviation of predicted states
        for j in range(grampc.opt.Nhor):
            points = np.reshape(grampc.rws.x[:, j], (Nx, point_transformation.number_of_points()), order='F')
            P = point_transformation.covariance(points, points)
            pred_mean[j, :] = point_transformation.mean(points)
            pred_std[j, :] = np.sqrt(np.diag(P))
        
        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t + grampc.param.dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.param.p0, grampc.param))
        next_state = sol.y[:Nx, -1]

        # update distribution of initial state
        initial_state.set_mean_and_covariance(next_state, covariance)
        stochastic_problem.compute_x0_and_p0(initial_state)
        # set current time and state
        grampc.set_param({"t0": t + grampc.param.dt, "x0": stochastic_problem.x0})
        
        # update plots
        if i % plotSteps == 0:
            ph_vec.set_data(vec.mean[:, 0], vec.mean[:, 1])
            for j in range(point_transformation.number_of_points()):
                ph_pred[j].set_data(grampc.rws.x[Nx*j, :], grampc.rws.x[Nx*j+1, :])
            for j in range(Nx):
                ph_vec_mean[j].set_data(vec.t, vec.mean[:, j])
                ph_vec_std_up[j].set_data(vec.t, vec.mean[:, j] + zalpha * vec.std[:, j])
                ph_vec_std_down[j].set_data(vec.t, vec.mean[:, j] - zalpha * vec.std[:, j])
                ph_pred_mean[j].set_data(t + grampc.rws.t, pred_mean[:, j])
                ph_pred_std_up[j].set_data(t + grampc.rws.t, pred_mean[:, j] + zalpha * pred_std[:, j])
                ph_pred_std_down[j].set_data(t + grampc.rws.t, pred_mean[:, j] - zalpha * pred_std[:, j])
            plt.pause(grampc.param.dt)

print(f'Computation time: mean {np.mean(vec.CPUtime):.3f} ms and max {np.max(vec.CPUtime):.3f} ms')

plt.ioff()
plt.show()