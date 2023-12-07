directory = '../../bin/';

% vec.J = load([directory 'Jvec.txt']);
vec.t = load([directory 'tvec.txt']);
vec.x = load([directory 'xvec.txt']);
vec.adj = load([directory 'adjvec.txt']);
vec.u = load([directory 'uvec.txt']);

T = vec.t(end);

vec.mean = vec.x(:, 1:6);
if(size(vec.x, 2) > 6)
    vec.cov = vec.x(:, 7:42);
    vec.std = sqrt(vec.cov(:, [1 8 15 22 29 36]));
else
    vec.std = zeros(size(vec.x, 1), 6);
end

figure(2);
subplot(2, 2, 1);
plot(vec.t, vec.mean(:, 4));
hold on
plot(vec.t, [vec.mean(:, 4) + vec.std(:, 4), vec.mean(:, 4) - vec.std(:, 4)], 'Color', [0.5 0.5 0.5 0.5]);
hold off
xlabel('t');
ylabel('E[x_1]');

subplot(2, 2, 2);
plot(vec.t, vec.mean(:, 5));
hold on
plot(vec.t, [vec.mean(:, 5) + vec.std(:, 5), vec.mean(:, 5) - vec.std(:, 5)], 'Color', [0.5 0.5 0.5 0.5]);
% plot([vec.t(1), vec.t(end)], [5, 5], 'k--');
hold off
xlabel('t');
ylabel('E[x_2]');

subplot(2, 2, 3);
plot(vec.t, vec.u);
xlabel('t');
ylabel('u');

subplot(2, 2, 4);
plot(vec.mean(:, 4), vec.mean(:, 5));
hold on
plot(vec.mean(:, 4) + vec.std(:, 4), vec.mean(:, 5) + vec.std(:, 5), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 4) + vec.std(:, 4), vec.mean(:, 5) - vec.std(:, 5), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 4) - vec.std(:, 4), vec.mean(:, 5) + vec.std(:, 5), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 4) - vec.std(:, 4), vec.mean(:, 5) - vec.std(:, 5), 'Color', [0.5 0.5 0.5 0.5]);
hold off
xlabel('x_1');
ylabel('x_2');
xlim([0 15*T]);
ylim([0 0.6]);

