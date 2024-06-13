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

% plot states
figure(1);
subplot(2, 1, 1);
plot(vec.t, vec.x);
xlabel('Time');
ylabel('States');
title('States');

% plot controls
subplot(2, 1, 2);
plot(vec.t(1:end-1), vec.u);
xlabel('Time');
ylabel('Controls');
title('Controls');


figure(2)
steps = 8;
plot(vec.x(1:steps:end,1), vec.x(1:steps:end,2),'.-','MarkerSize',10);
hold on;
plot([-50;50], [6 0;6 0],'k--');
pos = [18 -1.5 4 4];
rectangle('Position',pos,'Curvature',[1 1])
hold off, 
axis equal
ylim([-10,10]),
xlim([-2,vec.x(end,1)+2])
xlabel('x_1'),ylabel('x_2')


