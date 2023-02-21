clear all; %#ok<CLALL>
close all;
clc;

%% Problem setup

L = 100; % signal length
K = 2;  % number of different signals to estimate (heterogeneity)

% Ground truth signals
% x_true = transpose(linspace(-1,1,L));
x_true = transpose(sin(linspace(0,2*pi,L)));
% x_true(11:20) = transpose(sin(linspace(0,4*pi,10)));
% x_true(end-9:end) = transpose(sin(linspace(0,4*pi,10)));
x_true = normalize(x_true);

% Ground truth mixing probabilities
p_true = rand(K, 1);
p_true = max(.2*(1/K), p_true);
p_true = p_true / sum(p_true);
% Number of measurements for each class
M = 1e5;
Ms = round(p_true*M);

% Noise level
sigma = 0.1;
shift = -1; % shift set to -1 to enable cyclic shifts
% Generate the data
data = generate_observations_het(x_true, Ms, sigma, shift);


%% Optimization

opts = struct();
opts.maxiter = 200;
opts.tolgradnorm = 1e-7;
opts.tolcost = 1e-18;

% initial point
X0 = zeros(L,1);
p0 = ones(K, 1) / K;

[x_est, p_est, problem] = MRA_het_mixed_invariants_free_p(data, sigma, K, [], [], opts);

%% Evaluate quality of recovery, up to permutations and shifts.

[ind, x_est, ~, perm] = align_to_reference_het(x_est, x_true);
p_est = p_est(perm);

rel_error_X = norm(x_est - x_true) / norm(x_true);
tv_error_p = norm(p_est - p_true, 1) / 2;


%% Plot
d1 = floor(sqrt(K));
d2 = ceil(K/d1);
for k = 1 : K
    subplot(d1, d2, k);
    plot(1:L, x_true(:, k), '.-', 1:L, x_est(:, k), 'o-');
    legend('True signal', 'Estimate');
    title(sprintf('True weight: %.3g; estimated: %.3g, Relative Error: %.3g, Lag: %i', p_true(k), p_est(k), rel_error_X(k), ind(k)));
end
