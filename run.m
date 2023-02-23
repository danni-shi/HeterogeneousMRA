function run(L, K, sigma, max_shift, signal, save_data)
% run the entire data generation, signal recovery and recovery evaluation 
% process with different experiment settings.
    
    % Default values
    arguments
        L double {mustBePositive} = 50;
        K double {mustBePositive} = 2;
        sigma double {mustBeNonnegative} = 1;
        max_shift double {mustBeLessThan(max_shift, 1)} = 0.1;
        signal string = 'gaussian';
        save_data logical = true;
    end

    % Ground truth signals
    if signal == 'gaussian'
        x_true = randn(L, K);
        x_true = normalize(x_true);
    end
    
    % Sine waves
    % for k = 1:K
    %     x_true(:,k) = transpose(sin(linspace(0,k*pi,L)));
    % end
    
    % Ground truth mixing probabilities
    p_true = rand(K, 1);
    p_true = max(.2*(1/K), p_true);
    p_true = p_true / sum(p_true);
    
    % Number of measurements for each class
    M = 1e5;
    Ms = round(p_true*M);
    
    % Generate the data
    [data, shifts, classes] = generate_observations_het(x_true, Ms, sigma, max_shift);
    
    
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
    
    [x_est, E, perm] = align_to_reference_het(x_est, x_true);
    p_est = p_est(perm);
    rel_error_X = zeros(K,1);
    tv_error_p = norm(p_est - p_true, 1) / 2;
    
    for k = 1 : K
        rel_error_X(k) = norm(x_est(:, k) - x_true(:, k)) / norm(x_true(:, k));
    end
    
    if save_data
        save(sprintf('data/observations_noise%.3g_shift%.3g_class%i.mat', sigma, max_shift, K), 'data', 'shifts', 'classes');
        save(sprintf('data/results_noise%.3g_shift%.3g_class%i.mat', sigma, max_shift, K), 'x_est', 'p_est', 'rel_error_X', 'tv_error_p')
    else              
        %% Plot
        d1 = floor(sqrt(K));
        d2 = ceil(K/d1);
        for k = 1 : K
            subplot(d1, d2, k);
            plot(1:L, x_true(:, k), '.-', 1:L, x_est(:, k), 'o-');
            legend('True signal', 'Estimate');
            title(sprintf('True weight: %.3g; estimated: %.3g, Relative Error: %.3g, Lag: %i', p_true(k), p_est(k), rel_error_X(k)));
        end
    end
    

