%% Blatt 6, Aufgabe 1(iii) - Felix Benning und Nina Siebke

%initialise given parameters
T = 2;
x_0 = 1;
mu = 2;
lambda = 0.7;
n = 10^3;

time = 0:T/n:T;
num_sims = 3;

elapsed = zeros(1,num_sims);
simulatons = zeros(num_sims, length(time));
for idx = 1:num_sims
    tic
    simulatons(idx,:) = ornsteinUhlenbeck(T, x_0, mu, lambda, n);
    elapsed(idx) = toc;
end

plot(time,simulatons);


%% monte carlo

% theoretical values
theoretical_mean = x_0.*exp(-lambda.*time)+mu*(1-exp(-lambda.*time));
theoretical_var = 1/(2*lambda)*(1-exp(-2*lambda.*time));

% monte carlo values
simulation_mean = 1/num_sims.*sum(simulatons,1);
simulaton_variance = 1/num_sims .* sum(simulatons.*simulatons,1);
Fehler = mean(abs(theoretical_mean-simulation_mean));
FehlerVar = mean(abs(simulaton_variance-theoretical_var)); 


