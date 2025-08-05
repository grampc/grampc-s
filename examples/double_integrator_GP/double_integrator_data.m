% This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
%
% GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
%
% Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
% All rights reserved.
%
% GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

clear all; clc; close all;

%% Generate input data
x2 = -2 : 0.5 : 2;
u = -1 : 0.4 : 1;
[X2, U] = meshgrid(x2, u);

inputData = [X2(:)'; U(:)'];
numDatapoints = size(inputData, 2);

%% Generate output data with noise
dynamic_1 = @(dataPoint) [dataPoint(1,:)];
dynamic_2 = @(dataPoint) [dataPoint(2,:)];

noiseStandardDeviation = 1e-1;
outputData_1 = dynamic_1(inputData) + noiseStandardDeviation * randn(1, numDatapoints);
outputData_2 = dynamic_2(inputData) + noiseStandardDeviation * randn(1, numDatapoints);

%% Write data to file
addpath('../../matlab/')

writeGaussianProcessData('../../bin/GP1.txt', inputData, outputData_1, noiseStandardDeviation);
writeGaussianProcessData('../../bin/GP2.txt', inputData, outputData_2, noiseStandardDeviation);