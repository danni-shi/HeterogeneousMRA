clear all; %#ok<CLALL>
close all;
clc;

%% Get data
% daily returns from close to close
data_table = readtable('../data/pvCLCL_clean.csv'); 
tickers = data_table(:,1);
data_table(:,1) = [];
period_length = 50;
period_retrain = 10; % retrain every 10 trading days

K_range = 1:3;
sigma_range = 0.2:0.2:2.0;

%% Optimization
opts = struct();
opts.maxiter = 200;
opts.tolgradnorm = 1e-7;
opts.tolcost = 1e-18;
opts.verbosity = 0;

nextrainits = 2;


Nk = length(K_range);
Ns = length(sigma_range);

starting = 506;
ending = 1180;

for start_index = starting:period_retrain:ending
    tic;
    end_index = start_index + period_length - 1;
    data = transpose(table2array(data_table(:,start_index:end_index))); % L by N
    data = normalize(data,1); % normalize each column
    % data = data./std(data,0,2); % instead of normalize we divide by row std
    x_est_results = cell(Nk*Ns,1);
    p_est_results = cell(Nk*Ns,1);
    
    parfor i = 1:Nk*Ns
        [k, s] = ind2sub([Nk, Ns], i);
        K = K_range(k);
        sigma = sigma_range(s);
        fprintf('K = %i, sigma = %.3g \n', K, sigma)
        [x_est, p_est, problem] = MRA_het_mixed_invariants_free_p(data, sigma, K, [], [], opts, [], nextrainits);
        x_est_results{i} = x_est;
        p_est_results{i} = p_est;
    end
    % saving results using a for loop
    folder_name = sprintf('../data/pvCLCL_results/start%i_end%i',start_index,end_index);
        if exist(folder_name, 'dir')
            rmdir(folder_name,'s');
        end
        mkdir(folder_name);
    
    for i = 1:Nk*Ns
        [k, s] = ind2sub([Nk, Ns], i);
        K = K_range(k);
        sigma = sigma_range(s);
        file_specs = sprintf('noise%.3g_class%i',sigma, K);
        x_est = x_est_results{i};
        p_est = p_est_results{i};
        save(sprintf('%s/results_%s.mat', folder_name, file_specs), 'x_est','p_est');
    end
    toc;
end
% start_index = find(categorical(data_table.Properties.VariableNames) == {start_date});
% tickers
% names = data_table(:,1);
% names = join(erase(string(names{:, :}), "'"), '', 2);

% tickers = {'XLF','XLB','XLK','XLV','XLI','XLU','XLY','XLP','XLE'};
% index = find(ismember(names, tickers));
% M = length(index);
% d1 = floor(sqrt(M));
% d2 = ceil(M/d1);
% 
% for m = 1:M
% %     n = index(m);
%     n = m;
%     x = data(:,n);
%     [ind, x, ~, perm] = align_to_reference_het(x, x_est);
%     rel_error_X = norm(x_est - x) / norm(x_est);
%     subplot(d1, d2, m);
%     plot(1:L, x, '.-', 1:L, x_est, 'o-');
%     legend('Observation', 'Estimate of Signal');
%     title(sprintf('Ticker: %s ; Relative Error: %.3g, Lag: %i', names(n), rel_error_X, ind));
% end
