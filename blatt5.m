%% special variable definitions

dim=1;
sigma =0.3;
T = 2;
s = 2.7;
r = 0.05;
alpha = 1/2;

epsilons = (ones(1,6)*2).^-(2:7);
len_eps = length(epsilons);
conf_niveau = 0.05;

functional = @(t,bm) functional_G(t,black_scholes(t,bm,sigma,r,s),r);

%% Std Monte Carlo

est = zeros(1,len_eps);
variance = zeros(1,len_eps);
conf_interval = zeros(2,len_eps);
times = zeros(1, len_eps);

rv_generator = @(t) functional(t, brownian_motion(dim,t));

for idx = 1:len_eps
    tic
    eps = epsilons(idx);
    n = ceil(eps^-2);
    [est(idx), variance(idx), conf_interval(:,idx)] = monte_carlo(rv_generator, n, alpha, T, conf_niveau);
    times(idx) = toc;
end

solution_array = cat(2,epsilons', est', variance', conf_interval', times');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'epsilon', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'standard_5.csv')
loglog(times, conf_interval(2,:) - conf_interval(1,:), 'o-')

%% Romberg Monte Carlo
beta = 1/2;

est = zeros(1,len_eps);
variance = zeros(1,len_eps);
conf_interval = zeros(2,len_eps);
times = zeros(1, len_eps);

for idx = 1:length(epsilons)
    tic
    eps = epsilons(idx);
    n = ceil((sqrt(2/3)*eps)^-1/alpha);
    [est(idx), variance(idx), conf_interval(:,idx)] = romberg_monte_carlo(functional, n, alpha,beta, T, dim, conf_niveau);
    times(idx) = toc;
end

solution_array = cat(2,epsilons',est', variance', conf_interval', times');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'epsilon', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'romberg_5.csv')
loglog(times, conf_interval(2,:) - conf_interval(1,:), 'o-')

%% Multilevel Monte Carlo

est = zeros(1,len_eps);
variance = zeros(1,len_eps);
conf_interval = zeros(2,len_eps);
times = zeros(1, len_eps);

for idx = 1:length(epsilons)
    tic
    eps = epsilons(idx);
    times(idx) = toc;
end

solution_array = cat(2,epsilons',est', variance', conf_interval', times');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'epsilon', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'multilevel_5.csv')
loglog(times, conf_interval(2,:) - conf_interval(1,:), 'o-')

%% general function definitions

function payoff = functional_G(time_series, price_series, interest)
    % 1xT dimensional time_series and price_series
    % :return payoff
    expiration = time_series(end);
    payoff = exp(-interest*expiration)*(price_series(expiration) - min(price_series));
end

function [estimator, variance, conf_interval] = multilevel_monte_carlo(functional)
    
end