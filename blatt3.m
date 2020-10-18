%% special variable definitions

dim=2;
sigma =[
    0.3 0.12
    0.12 0.2
];
T = 1;
s = [3.3; 2.1];
r = 0.06;

%% a)

alpha = 0.5;
beta = 1/2;
n_vector = [10, 10^2, 10^3, 10^4, 10^5];
niveau = 0.05;

estimator = zeros(1,length(n_vector));
variance = zeros(1,length(n_vector));
conf_interval = zeros(2,length(n_vector));
times = zeros(1, length(n_vector));

for i=1:length(n_vector)
    tic
    functional = @(t, bm) G(black_scholes(t, bm, sigma, r, s),T,r);
    [estimator(i), variance(i), conf_interval(:,i)] = romberg_monte_carlo(functional, n_vector(i), alpha, beta, T, dim, niveau);
    times(i) = toc;
end

solution_array = cat(2,n_vector',estimator', variance', conf_interval', times');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'N', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'romberg.csv')
plot(times, estimator, 'o-')

%% general function definitions

function [estimator,variance,conf_interval] = romberg_monte_carlo(functional, n, alpha, beta, T, dim, niveau)
    % N: number of samples
    % niveau: niveau for confidence_interval
    % rv_generator: function which generates 1-dim random variables
    % :return: monte carlo estimator, estimated variance and confidence
    %          interval
    m = ceil(n^beta);
    M = ceil(n^(2*alpha));
    N = ceil(n^(2*alpha-beta));
    times_m = 0:1/m:T;
    
    p_M = zeros(M,1);
    for idx = 1:M
        bm = brownian_motion(dim, times_m);
        p_M(idx) = functional(times_m, bm);
    end
    
    
    times_diff = cell(2,1);
    times_diff{1} = 0:1/n:T;
    times_diff{2} = times_m;
    
    p_diff = zeros(N,1);
    for idx = 1:N
        bm = brownian_motion(dim, times_diff);
        p_diff(idx) = functional(times_diff{1}, bm{1}) - functional(times_diff{2}, bm{2});
    end
    
    estimator = mean(p_M) + mean(p_diff);
    variance = var(p_M)/M + var(p_diff)/N;
    conf_interval = confidence_interval(estimator, variance, niveau);
end
   