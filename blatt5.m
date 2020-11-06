%% special variable definitions

dim=1;
sigma =0.3;
T = 2;
s = 2.7;
r = 0.05;
alpha = 1/2;

epsilons = arrayfun(@(x) 2^-x, 2:7);
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
M = 4;
c_var = 1/2;

est = zeros(1,len_eps);
variance = zeros(1,len_eps);
conf_interval = zeros(2,len_eps);
times = zeros(1, len_eps);

for idx = 1:length(epsilons)
    tic
    eps = epsilons(idx);
    L = ceil(log(1/eps)/log(2));
    N = arrayfun(@(l) ceil(2/(eps^2)*(L+1)*c_var*T*M^-(l-1)) ,1:(L+1));
    [est(idx), variance(idx), conf_interval(:,idx)] = multilevel_monte_carlo(functional, dim, T, M, L, N, conf_niveau); 
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
    payoff = exp(-interest*expiration)*(price_series(end) - min(price_series));
end

function [estimate, variance, conf_interval] = multilevel_monte_carlo(functional, dim, T, M, L, N, niveau)
    % first level
    times = cell(2,1);
    % n = 1; %M^0
    times{1} = [0,T]; % 0:T/n:T
    
    p_N = cell(L+1,1);
    p_N{1} = zeros(N(1),1);
    for idx = 1:N(1)
        bm = brownian_motion(dim, times{1});
        p_N{1}(idx) = functional(times{1}, bm); 
    end
    
    % remaining levels
    for level = 1:L
        times{2} = times{1};
        n = M^level;
        times{1} = 0:T/n:T;
        
        p_level = zeros(N(level+1),1);
        for idx = 1:N(level+1)
            bm = brownian_motion(dim, times);
            p_level(idx) = functional(times{1}, bm{1}) - functional(times{2}, bm{2});
        end
        p_N{level+1} = p_level;
    end
    
    % calculate mean, variance, confidence intervals
    estimate = 0;
    variance = 0;
    for level = 1:L+1
        estimate = estimate + mean(p_N{level});
        variance = variance + var(p_N{level})/N(level);
    end
    conf_interval = confidence_interval(estimate, variance, niveau);
end