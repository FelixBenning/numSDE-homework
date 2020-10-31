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
    {'n', 'mean', 'variance', '95% interv lower', '95% interv upper', 'time'}...
);
writetable(soultion_table, 'romberg.csv')
plot(times, estimator, 'o-')



