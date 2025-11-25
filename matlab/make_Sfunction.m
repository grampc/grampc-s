% This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
%
% GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
%
% Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
% All rights reserved.
%
% GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

function make_Sfunction(path_to_grampc_s, probfct_inc, probfct_src, s_function)
%MAKE_SFUNCTION Compile a S-Function for GRAMPC-S
%   Inputs: 
%      path_to_grampc_s : Path to GRAMPC-S
%      probfct_inc      : Source file of the problem description including the path to the file
%      probfct_src      : Header file of the problem description including the path to the file
%      s_function       : S-Function file including the path to the file

%% Paths
% Paths to C files
cFiles = [fullfile(path_to_grampc_s, "libs", "grampc", "src", ["discrete.c", "finite_diff.c", "grampc_alloc.c", "grampc_fixedsize.c", "grampc_init.c", ...
    "grampc_mess.c", "grampc_run.c", "grampc_setopt.c", "grampc_setparam.c", "grampc_util.c", "grampc_erk.c", "rodas.c", "ruku45.c", "simpson.c", ...
    "trapezoidal.c"])];

objFilesC = [fullfile("build", ["discrete.obj", "finite_diff.obj", "grampc_alloc.obj", "grampc_fixedsize.obj", "grampc_init.obj", ...
    "grampc_mess.obj", "grampc_run.obj", "grampc_setopt.obj", "grampc_setparam.obj", "grampc_util.obj", "grampc_erk.obj", "rodas.obj", "ruku45.obj", "simpson.obj", ...
    "trapezoidal.obj"])];

% Paths to C++ files
cppFiles = [fullfile(path_to_grampc_s, "src", "constraint_approx", ["chebyshev_constraint_approximation.cpp", "Gaussian_constraint_approximation.cpp", "symmetric_constraint_approximation.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "distribution", ["distribution.cpp", "Gaussian_distribution.cpp", "multivariate_uncorrelated_distribution.cpp", "univariate_uniform_distribution.cpp", ...
    "univariate_beta_distribution.cpp", "univariate_chi_squared_distribution.cpp", "univariate_extreme_value_distribution.cpp", ...
    "univariate_f_distribution.cpp", "univariate_gamma_distribution.cpp", "univariate_lognormal_distribution.cpp", "univariate_piecewise_constant_distribution.cpp", "univariate_student-t_distribution.cpp", ...
    "univariate_weibull_distribution.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "gaussian_process", ["gaussian_process.cpp", "kernel_product.cpp", "kernel_sum.cpp", "locally_periodic_kernel.cpp", "Matern32_kernel.cpp", "Matern52_kernel.cpp", ...
    "periodic_kernel.cpp", "squared_exponential_kernel.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "grampc_interface", "grampc_interface.cpp"), ...
    fullfile(path_to_grampc_s, "src", "point_transformation", ["composed_quadrature.cpp", "Hermite_quadrature.cpp", "Legendre_quadrature.cpp", "monte_carlo.cpp", "PCE_transformation.cpp", ...
    "stirling_interpolation_first_order.cpp", "stirling_interpolation_second_order.cpp", "unscented_transformation.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "polynomial", ["Hermite_polynomial_generator.cpp", "Legendre_polynomial_generator.cpp", "multivariate_polynomial.cpp", "polynomial.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "problem_description", ["sigma_point_problem_description.cpp", "monte_carlo_problem_description.cpp", "resampling_GP_problem_description.cpp", ...
    "resampling_problem_description.cpp", "taylor_problem_description.cpp", "problem_description.cpp"]), ...
    fullfile(path_to_grampc_s, "src", "simulator", "simulator.cpp"), ...
    fullfile(path_to_grampc_s, "src", "util", "grampc_s_util.cpp")];

objFilesCpp = [fullfile("build", ["chebyshev_constraint_approximation.obj", "Gaussian_constraint_approximation.obj", "symmetric_constraint_approximation.obj", ...
    "distribution.obj", "Gaussian_distribution.obj", "multivariate_uncorrelated_distribution.obj", "univariate_uniform_distribution.obj", ...
    "univariate_beta_distribution.obj", "univariate_chi_squared_distribution.obj", "univariate_extreme_value_distribution.obj", ...
    "univariate_f_distribution.obj", "univariate_gamma_distribution.obj", "univariate_lognormal_distribution.obj", "univariate_piecewise_constant_distribution.obj", "univariate_student-t_distribution.obj", ...
    "univariate_weibull_distribution.obj", "gaussian_process.obj", "kernel_product.obj", "kernel_sum.obj", "locally_periodic_kernel.obj", "Matern32_kernel.obj", "Matern52_kernel.obj", ...
    "periodic_kernel.obj", "squared_exponential_kernel.obj", "grampc_interface.obj", "composed_quadrature.obj", "Hermite_quadrature.obj", "Legendre_quadrature.obj", ...
    "monte_carlo.obj", "PCE_transformation.obj", "stirling_interpolation_first_order.obj", "stirling_interpolation_second_order.obj", "unscented_transformation.obj",  ...
    "Hermite_polynomial_generator.obj", "Legendre_polynomial_generator.obj", "multivariate_polynomial.obj", "polynomial.obj", ...
    "sigma_point_problem_description.obj", "monte_carlo_problem_description.obj", "resampling_GP_problem_description.obj", ...
    "resampling_problem_description.obj", "taylor_problem_description.obj", "problem_description.obj", "simulator.obj", "grampc_s_util.obj"])];

% Include folders
[probfct_dir,~,~] = fileparts(probfct_inc);
Ipath1 = ['-I' fullfile(path_to_grampc_s, 'libs', 'grampc', 'include') ' '];
Ipath2 = ['-I' probfct_dir ' '];
Ipath3 = ['-I' fullfile(path_to_grampc_s, 'include') ' '];
if ispc
    Ipath4 = ['-I' '"C:\Program Files (x86)\Eigen3\include\eigen3"' ' '];
else
    Ipath4 = ['-I' '/usr/local/include/eigen3' ' '];
end
%% Compile c and c++ files

% Generate string with all C files
Objc = '';
for i = 1:length(cFiles)
    Objc = [Objc convertStringsToChars(cFiles(i))  ' '];
end

% Generate string with all C obj-files
Objc_compiled = '';
for i = 1:length(cFiles)
    Objc_compiled = [Objc_compiled convertStringsToChars(objFilesC(i))  ' '];
end

% Generate string with all C++ files
Objcpp = '';
for i = 1:length(cppFiles)
    Objcpp = [Objcpp convertStringsToChars(cppFiles(i))  ' '];
end

% Generate string with all C obj-files
Objcpp_compiled = '';
for i = 1:length(cppFiles)
    Objcpp_compiled = [Objcpp_compiled convertStringsToChars(objFilesCpp(i))  ' '];
end

% Path to S-function
[~, s_function_name, ~] = fileparts(s_function);

% Build C and C++ files
disp('Building S-function ...');
if ~isfolder('build')
    eval(['mex -v -c -DMEXCOMPILE -outdir build ' Objc ' ' Ipath1]);
    eval(['mex -v -c -outdir build ' char("CXXFLAGS='$CXXFLAGS -std=c++11'") ' -DMEXCOMPILE ' Objcpp Objc_compiled ' ' Ipath1 Ipath2 Ipath3 Ipath4]);
end
eval(['mex -v ' char("CXXFLAGS='$CXXFLAGS -std=c++11'") ' -output ' s_function_name ' -DMEXCOMPILE ' s_function ' ' Objc_compiled ' ' Objcpp_compiled ' ' probfct_src ' ' Ipath1 Ipath2 Ipath3 Ipath4]);
disp('S-function is compiled successfully.');
end