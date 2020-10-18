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
        prices = black_scholes(dim, T, sigma, r, s);
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
    times = 0:1/n_vector(i):T;
    rv_generator = @()G(black_scholes(dim, times, sigma, r, s),T,r);
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
   


