# GRAMPC-S
GRAMPC-S is a framework for nonlinear stochastic model predictive control based on [GRAMPC](https://github.com/grampc/grampc). The implementation in C++ and the fast solver allow sampling times in the (sub)millisecond range. A manual of GRAMPC-S is provided in the folder [documentation](documentation/manual.pdf).

More details about the algorithm and its performance can be found in the corresponding article published in Optimization and Engeneering. The article is available online under open access at: https://doi.org/10.1007/s11081-025-10006-z. Please cite the paper when you are using results obtained with GRAMPC-S.

## Features
The most important features are listed below. For a more detailed description of the features, please refer to the [documentation](documentation/manual.pdf).
- Approximation of a stochastic optimal control problem by a deterministic optimal control problem that can be solved by [GRAMPC](https://github.com/grampc/grampc),
- Several approximation methods for predicting stochastic moments of random variables:
  - First-order Taylor linearization,
  - Stochastic sampling,
  - Transformations between points and stochastic moments, such as the unscented transformation, multiple Gaussian quadrature rules, and Stirling's interpolation of first and second order,
  - Polynomial chaos expansion.
- Approximation of probabilistic inequality constraints by deterministic inequality constraints,
- Random sampling from a large number of probability density functions,
- Consideration of additive Gaussian processes in the system function for modeling unknown parts of the dynamics using measured data points.

## Dependencies
GRAMPC-S depends on the [Eigen library](https://eigen.tuxfamily.org) in version 3.4.90 or higher (tested with version 3.4.90). If Eigen is not installed on your system, clone the repository, create the folder *build*, and install Eigen. For **Linux**, execute the following commands:
```
git clone https://gitlab.com/libeigen/eigen.git
mkdir eigen/build
cd eigen/build
cmake -DCMAKE_BUILD_TYPE=Release ..
sudo make install
```
For **Windows** with MinGW compiler, execute the following commands:
```
git clone https://gitlab.com/libeigen/eigen.git
mkdir eigen\build
cd eigen\build
cmake -DCMAKE_BUILD_TYPE=Release .. -G "MinGW Makefiles"
cmake --install .
```

## Installation
Clone the repository and all submodules:
```
git clone --recurse-submodules https://github.com/grampc/grampc-s.git
```
To build GRAMPC-S using Linux, run the CMakeLists.txt and call make afterwards:
```
cmake -DCMAKE_BUILD_TYPE=Release CMakeLists.txt
make
```
To build GRAMPC-S using Windows, the editor Visual Studio Code with the [C++ extension](https://code.visualstudio.com/docs/languages/cpp) is recommended. Open the folder "grampc-s" in Visual Studio Code, click the button "Build", and select a compiler from the list.

## License
GRAMPC-S is licensed under BSD 3-clause license, see LICENSE.txt.