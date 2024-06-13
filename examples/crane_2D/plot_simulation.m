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

%% Additional code
figure(2),clf;
ph = subplot(1,1,1);
axis equal
steps = 150;
x_temp = vec.x(floor(1:numel(vec.x(:,1))/steps:numel(vec.x(:,1))),:)';
yy=-2.5:0.01:2.5;
yyy = -0.2*yy.^2-1.25;
plot(yy,yyy,'Parent',ph)
hold on
for i = 1:numel(x_temp(1,:))
    %         tic;
    draw_CRANE_2D(x_temp(:,i), 0.3+0.7/(steps+2)*i,ph)
    drawnow
    pause(0.01);
end
axis(ph,'equal')

function draw_CRANE_2D(x, color, ph)
temp_color = [1-color, 1-color, 1-color];
rectangle('Position',[x(1)-0.25 -0.1 0.5 0.2],'Curvature',0.2,'EdgeColor',temp_color,'Parent',ph)
rectangle('Position',[x(1)-0.22 -0.12 0.08 0.08],'Curvature',[1 1],'EdgeColor',temp_color,'Parent',ph)
rectangle('Position',[x(1)+0.14 -0.12 0.08 0.08],'Curvature',[1 1],'EdgeColor',temp_color,'Parent',ph)
plot([x(1) x(1)+(x(3)-0.08)*sin(x(5))],[0, -(x(3)-0.08)*cos(x(5))],'Color',temp_color,'Parent',ph)
rectangle('Position',[x(1)+x(3)*sin(x(5))-0.04 -x(3)*cos(x(5))+0.00 0.08 0.08],'Curvature',[1 1],'EdgeColor',temp_color,'Parent',ph)
end
