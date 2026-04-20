# Installation

Install Python interface for GRAMPC by following the steps in the [documentation](https://grampc.github.io/grampc/install.html#installation-of-grampc-for-use-in-python).

Install Python interface for GRAMPC-S by executing ```pip install .``` in the root folder grampc-s.

Install mobile robot interface by executing ```pip install .``` in the subfolder grampc-s/python/examples/mobile_robot.

# Running

Run ```python mobile_robot_with_uncertain_initial_state.py``` for simulating the mobile robot with uncertainty in the initial state.

Run `python mobile_robot_with_process_noise.py` for simulating the mobile robot with noise in the dynamics.

Run `python mobile_robot_with_uncertain_parameters.py` for simulating the mobile robot with uncertainty in the parameters that scale the linear and angular velocity.

# Troubleshooting

If compilation fails, check that you have installed the Python interfaces of GRAMPC and GRAMPC-S using the same Python version and (virtual) environment.

If the simulation results are not displayed, install PyQt by `pip install PyQt6`.
