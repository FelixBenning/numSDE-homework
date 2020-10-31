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
        prices = black_scholes(T, brownian_motion(dim,T), sigma, r, s);
        payoffs(k) = basket_call(prices, weights, strike, T, r);
    end
    avg(n) = mean(payoffs);
    variance(n) = var(payoffs)/N;
    conf_interval(:,n) = confidence_interval(avg(n), variance(n), niveau);
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
times = zeros(1, length(n_vector));

for i=1:length(n_vector)
    tic
    rv_generator = @(t) G(black_scholes(t, brownian_motion(dim,t), sigma, r, s),T,r);
    [estimator(i), variance(i), conf_interval(:,i)] = monte_carlo(rv_generator, n_vector(i), alpha, T, niveau);
    times(i) = toc;
end

solution_array = cat(2,n_vector',estimator', variance', conf_interval', times');
soultion_table = array2table(...
    solution_array,...
    'VariableNames',... 
    {'n', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'b_table.csv')
plot(times, estimator, 'o-')


