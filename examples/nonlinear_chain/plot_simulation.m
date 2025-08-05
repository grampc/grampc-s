% This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
%
% GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
%
% Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
% All rights reserved.
%
% GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
clear all; close all;

PLOT_MOVE_STEPS = 4;


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
N = (vec.Nx + 3) / 6;
subplot(5,1,4)
plot(vec.t(1:end-1), vec.u(:,1));
hold on;
plot(vec.t(1:end-1), vec.u(:,2));
plot(vec.t(1:end-1), vec.u(:,3));
title('controls');
subplot(5,1,5)
temp_vNorm = zeros(length(vec.t),1);
for ii = 1:length(vec.t)
    vVec = [];
    for i = 1:N-1
       vVec = [vVec vec.x(ii,(4+(i-1)*6):i*6)];
    end
    temp_vNorm(ii) = norm(vVec);
end
plot(vec.t, temp_vNorm);
title('||v||');
for ii = 1:PLOT_MOVE_STEPS:length(vec.t) - 2
    x1 = [0, vec.x(ii, 1:6:end)]';
    x2 = [0, vec.x(ii+1, 1:6:end)]';
    y1 = [0, vec.x(ii, 2:6:end)]';
    y2 = [0, vec.x(ii+1, 2:6:end)]';
    z1 = [0, vec.x(ii, 3:6:end)]';
    z2 = [0, vec.x(ii+1, 3:6:end)]';


    for jj = 0:2:9

        dx = x1 + jj * 0.1 * (x2 - x1);
        dy = y1 + jj * 0.1 * (y2 - y1);
        dz = z1 + jj * 0.1 * (z2 - z1);
        dux = vec.u(ii+1,1); 
        duy = vec.u(ii+1,2); 
        duz = vec.u(ii+1,3); 
        subplot(5,1,[1,2,3])
        hold off
        plot3(dx, dy, dz)
        grid on
        hold on
        plot3(dx, dy, dz, 'x', 'LineWidth',10)
        xlim([-1, N*0.6])
        ylim([-1, 1])
        zlim([-N*1.2 - 2, 4])
        
        title(['time = ', num2str(vec.t(ii))]);

        drawnow
    end
end
