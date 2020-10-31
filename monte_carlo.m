function [estimator,variance,conf_interval] = monte_carlo(rv_generator, n, alpha, T, niveau)
    % n, alpha: sampling parameters
    % niveau: niveau for confidence_interval
    % rv_generator: function which generates 1-dim random variables
    % :return: monte carlo estimator, estimated variance and confidence
    %          interval
    N = ceil(n^(2*alpha));
    times = 0:1/n:T;
    X_vector = arrayfun(@(variable) rv_generator(times),1:N);
    estimator = mean(X_vector);
    variance = var(X_vector)/N;
    conf_interval = confidence_interval(estimator, variance, niveau);
end
   