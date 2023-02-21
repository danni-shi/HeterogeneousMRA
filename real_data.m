clear all; %#ok<CLALL>
close all;
clc;

%% Get data
% daily returns from open to close
data_table = readtable('OPCL_20000103_20201231.csv'); 
% drop rows with missing data
data_table = data_table(~any(ismissing(data_table),2),:);
start_date = 'X20200102';
start_index = find(categorical(data_table.Properties.VariableNames) == {start_date});
% tickers
names = data_table(:,1);
names = join(erase(string(names{:, :}), "'"), '', 2);

% set the length of signal
L = 100;
sigma = 0;
K = 1;

data = transpose(table2array(data_table(:,start_index:start_index + L - 1)));
data = normalize(data,1);

%% Optimization
opts = struct();
opts.maxiter = 200;
opts.tolgradnorm = 1e-7;
opts.tolcost = 1e-18;
opts.verbosity = 1;

[x_est, p_est, problem] = MRA_het_mixed_invariants_free_p(data, sigma, K, [], [], opts);
p_est = p_est(perm);

tickers = {'XLF','XLB','XLK','XLV','XLI','XLU','XLY','XLP','XLE'};
index = find(ismember(names, tickers));
M = length(index);
d1 = floor(sqrt(M));
d2 = ceil(M/d1);

for m = 1:M
%     n = index(m);
    n = m;
    x = data(:,n);
    [ind, x, ~, perm] = align_to_reference_het(x, x_est);
    rel_error_X = norm(x_est - x) / norm(x_est);
    subplot(d1, d2, m);
    plot(1:L, x, '.-', 1:L, x_est, 'o-');
    legend('Observation', 'Estimate of Signal');
    title(sprintf('Ticker: %s ; Relative Error: %.3g, Lag: %i', names(n), rel_error_X, ind));
end
