% This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
%
% GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
%
% Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
% All rights reserved.
%
% GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

% directory with text files
directory = '../../bin/';

% load data
vec.t = load([directory 'tvec.txt']);
vec.x = load([directory 'xvec.txt']);
vec.u = load([directory 'uvec.txt']);
vec.dim = readtable([directory 'dim.txt'], "Delimiter", ",");

% number of states and controls
vec.Nx = vec.dim.Data(1);
vec.Nu = vec.dim.Data(2);
vec.Np = sqrt(size(vec.x, 2) - vec.Nx) - vec.Nx; 

% mean and standard deviation of the states
vec.mean = vec.x(:, 1 : vec.Nx);
vec.cov = vec.x(:, vec.Nx+1 : end);
vec.std = sqrt(vec.cov(:, 1 : (vec.Nx + vec.Np + 1) : end));

% plot states
figure(1);
subplot(2, 1, 1);
hold on;
for i = 1:vec.Nx
    line = plot(vec.t, vec.mean(:,i));
    patch([vec.t; flipud(vec.t)], [vec.mean(:,i) - vec.std(:,i); flipud(vec.mean(:,i) + vec.std(:,i))], line.Color,'FaceAlpha', 0.2, 'LineStyle', 'none');
end
hold off;
xlabel('Time');
ylabel("States");
title('Mean and standard deviation of states');

% plot controls
subplot(2, 1, 2);
plot(vec.t, vec.u);
xlabel('Time');
ylabel('Controls');
title('Controls');
