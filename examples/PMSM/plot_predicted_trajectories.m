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

% plot state trajectories
fig = figure('WindowState', 'maximized');
subplot(2,1,1);
colorOrder = get(gca, 'ColorOrder');
hold on
for i = 1 : vec.Nx
    plot(vec.t, vec.x(:, i:vec.Nx:end), 'Color', [colorOrder(mod(i-1, size(colorOrder, 1)) + 1, :), 1]);
end
hold off
xlabel('Time');
ylabel('States');
title('States')

% plot controls
subplot(2,1,2);
hold on
for i = 1:vec.Nu
    plot(vec.t, vec.u(:, i));
end
hold off
xlabel('Time');
ylabel('Controls');
title('Controls')

