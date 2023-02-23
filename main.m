clear all; %#ok<CLALL>
close all;
clc;

rng(1);
save_data = true;

%% Problem setup
L = 50; % signal length
% K = 1;  % number of different signals to estimate (heterogeneity)
signal = 'gaussian';
max_shift = 0.1; % shift set to -1 to enable cyclic shifts

%% generate data at different noise level and different number of classes
for K = 1:4
    for sigma = 0.1:0.1:2
        run(L, K, sigma, max_shift, signal, save_data);
    end
end
