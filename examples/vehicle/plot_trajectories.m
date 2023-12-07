clear all; close all; clc;
directory = '../../bin/';

% Plot trajectories of these states
plotState = 1:6;

% load data
vec.t = load([directory 'tvec.txt']);
vec.x = load([directory 'xvec.txt']);
vec.u = load([directory 'uvec.txt']);
vec.adj = load([directory 'adjvec.txt']);
vec.con = load([directory 'constr.txt']);
vec.dim = readtable([directory 'dim.txt'], "Delimiter", ",");

% prediction horizon
T = vec.t(end);

% dimensions
vec.Nx = vec.dim.Data(1);
vec.Nu = vec.dim.Data(2);
vec.Np = vec.dim.Data(3);
vec.Nh = vec.dim.Data(4);

% state mean
vec.mean = zeros(length(vec.t), vec.Nx);
for i = 1 : length(vec.t)
    for j = 1 : vec.Nx
        vec.mean(i,j) = mean(vec.x(i, j:vec.Nx:end));
    end
end

% plot states and mean
fig = figure('WindowState', 'maximized');
subplot(1, 3, 1);
colorOrder = get(gca, 'ColorOrder');
hold on
for i = 1 : length(plotState)
    plot(vec.t, vec.x(:, plotState(i):vec.Nx:end), 'Color', [colorOrder(mod(i-1, size(colorOrder, 1)) + 1, :), 0.2], 'HandleVisibility','off');
    plot(vec.t, vec.mean(:, plotState(i)), 'Color', colorOrder(mod(i-1, size(colorOrder, 1)) + 1, :), 'DisplayName', ['State ', num2str(plotState(i))]);
end
hold off
xlabel('Time');
ylabel('States');
title('States')
if length(plotState) > 1
    legend
end

% plot u
subplot(1, 3, 2);
set(gca, 'ColorOrderIndex', 1);
hold on
for i = 1:vec.Nu
    plot(vec.t, vec.u(:, i), 'Color', colorOrder(1,:), 'DisplayName', ['Control ', num2str(i)]);
end
hold off
xlabel('Time');
ylabel('Controls');
title('Controls')
if vec.Nu > 1
    legend
end

% plot constraints
subplot(1, 3, 3);
set(gca, 'ColorOrderIndex', 1);
hold on
plot([vec.t(1) vec.t(end)], [0 0], 'k--')
for i = 1 : vec.Nh
    plot(vec.t, vec.con(:, i), 'Color', [colorOrder(mod(i-1, size(colorOrder, 1)) + 1, :)]);
end
hold off
xlabel('Time');
ylabel('Constraints');
title('Constraints')

