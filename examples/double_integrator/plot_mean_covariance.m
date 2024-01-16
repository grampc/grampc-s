directory = '../../bin/';

% vec.J = load([directory 'Jvec.txt']);
vec.t = load([directory 'tvec.txt']);
vec.x = load([directory 'xvec.txt']);
vec.adj = load([directory 'adjvec.txt']);
vec.u = load([directory 'uvec.txt']);

T = vec.t(end);

vec.mean = vec.x(:, 1:2);
if(size(vec.x, 2) > 2)
    vec.cov = vec.x(:, 3:end);
    % vec.std = sqrt(vec.cov(:, [1 8 15 22 29 36]));
    vec.std = sqrt(vec.cov(:, [1 4]));
else
    vec.std = zeros(size(vec.x, 1), 6);
end

figure(6);
subplot(2, 2, 1);
plot(vec.t, vec.mean(:, 1));
hold on
plot(vec.t, [vec.mean(:, 1) + vec.std(:, 1), vec.mean(:, 1) - vec.std(:, 1)], 'Color', [0.5 0.5 0.5 0.5]);
hold off
xlabel('t');
ylabel('E[x_1]');

subplot(2, 2, 2);
plot(vec.t, vec.mean(:, 2));
hold on
plot(vec.t, [vec.mean(:, 2) + vec.std(:, 2), vec.mean(:, 2) - vec.std(:, 2)], 'Color', [0.5 0.5 0.5 0.5]);
% plot([vec.t(1), vec.t(end)], [5, 5], 'k--');
hold off
xlabel('t');
ylabel('E[x_2]');

subplot(2, 2, 3);
plot(vec.t, vec.u);
xlabel('t');
ylabel('u');

subplot(2, 2, 4);
plot(vec.mean(:, 1), vec.mean(:, 2));
hold on
plot(vec.mean(:, 1) + vec.std(:, 1), vec.mean(:, 2) + vec.std(:, 2), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 1) + vec.std(:, 1), vec.mean(:, 2) - vec.std(:, 2), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 1) - vec.std(:, 1), vec.mean(:, 2) + vec.std(:, 2), 'Color', [0.5 0.5 0.5 0.5]);
plot(vec.mean(:, 1) - vec.std(:, 1), vec.mean(:, 2) - vec.std(:, 2), 'Color', [0.5 0.5 0.5 0.5]);
hold off
xlabel('x_1');
ylabel('x_2');
xlim([0 2]);
ylim([-1 1]);

