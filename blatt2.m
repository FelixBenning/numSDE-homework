%% special variable definitions

dim=2;
sigma =[
    0.3 0.12
    0.12 0.2
];
T = 1;
s = [3.3; 2.1];
r = 0.06;

%% monte carlo
%% a)

simulations = [10^2, 10^3, 10^4, 10^5, 10^6];
strike = 2.5; % of the basket call
weights = [0.6, 0.4]; % of the basket call
niveau = 0.05; % of the confidence interval

avg = zeros(1,length(simulations));
variance = zeros(1,length(simulations));
conf_interval = zeros(2,length(simulations));
for n = 1:length(simulations)
    N=simulations(n);
    payoffs = zeros(1,N);
    for k = 1:N
        % N simulations stored in payoffs
        prices = gen_black_scholes(dim, T, sigma, r);
        payoffs(k) = basket_call(prices, weights, strike, T, r);
    end
    avg(n) = mean(payoffs);
    variance(n) = var(payoffs);
    conf_interval(:,n) = confidence_interval(avg(n), variance(n), niveau, N);
end

solution_array = cat(2,simulations',avg', variance', conf_interval');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'N', 'mean', 'variance', '95% interv lower', '95% interv upper'}...
);
writetable(soultion_table, 'a_table.csv')

%% b)
alpha = 0.5;
n_vector = [10, 10^2, 10^3, 10^4, 10^5];
niveau = 0.05;
dim = 2;
estimator = zeros(1,length(n_vector));
variance = zeros(1,length(n_vector));
conf_interval = zeros(2,length(n_vector));

for i=1:length(n_vector)
    N = ceil(n_vector(i)^(2*alpha));
    times = 1/n_vector(i)*[0:T];
    rv_generator = @()G(gen_black_scholes(dim, times, sigma, r),T,r);
   [estimator(i), variance(i), conf_interval(:,i)] = monte_carlo(N,rv_generator,niveau);
end

solution_array = cat(2,n_vector',estimator', variance', conf_interval');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'N', 'mean', 'variance', '95% interv lower', '95% interv upper'}...
);
writetable(soultion_table, 'b_table.csv')


%% general function definitions

function conf_interval = confidence_interval(avg, variance, niveau, N)
    conf_interval = avg+norminv([niveau/2, 1-niveau/2])*sqrt(variance/N);
end

function brownian_motion = generate_bm(dim, times)
    % dim: how many dimensions?
    % times: sorted row vector (1xT) of positive timestamps to sample the bm
    % :returns borwnian_motion: (dim x T) vector where every column i is the
    %   bm at time times(i)
    s = size(times);
    std_normal = randn(dim, s(2)); % s(2) = T (1:columsize, 2:rowsize)
    sqrt_diffs = sqrt(diff(cat(2,0, times))); % insert t_0 = 0
    brownian_motion = zeros(dim,s(2)); % preallocating the correct size
    for k = 1:s(2)
        brownian_motion(:,k) = std_normal(:,1:k) * transpose(sqrt_diffs(1,1:k));
    end
end

function price_series = black_scholes(times, brownian_motion, cov_matrix, interest)
    % times: sorted row vector (1xT) of positive timestamps for prices
    % brownian_motion: d x T brownian motion where column i matches times(i)
    % cov_matrix: a m x m covariance matrix to use
    % interest: the interest rate
    % :return price_series: a m x T price series is returned 
    vec = (interest-0.5*sum(cov_matrix.^2, 2));
    price_series = exp(vec * times + cov_matrix * brownian_motion);
end

function price_series = gen_black_scholes(dim, times, cov_matrix, interest)
    price_series = black_scholes(times, generate_bm(dim,times), cov_matrix, interest);
end

function payoff = basket_call(price_at_expiration, weights, strike, expiration, interest)
    % price_at_expiration: (dx1) price vector of the d products in the
    % basket
    % weights: (1xd) weights of the d products in the basket
    % strike: cost of the entire basket fixed at time 0
    % expiration: expiration time
    % interest: interest rate
    % :return payoff: payoff of the basket call
    payoff = exp(-interest*expiration)*weights*price_at_expiration - strike;
    payoff = max(payoff, 0);
end

function payoff = G(price_series, expiration, interest)
    % price_series: (2xT) dimensional price series
    % expiration: expiration time
    % interest: interest rate
    % :return payoff: payoff of this payoff function
    payoff = exp(-interest*expiration)*(max(price_series(1,:))-min(price_series(2,:)));
    payoff = max(payoff, 0);
end

function [estimator,variance,conf_interval] = monte_carlo(N, rv_generator, niveau)
    % N: number of samples
    % niveau: niveau for confidence_interval
    % rv_generator: function which generates 1-dim random variables
    % :return: monte carlo estimator, estimated variance and confidence
    %          interval
    X_vector = arrayfun(@(variable)rv_generator(),1:N);
    estimator = mean(X_vector);
    variance = var(X_vector);
    conf_interval = confidence_interval(estimator, variance, niveau, N);
end
   


