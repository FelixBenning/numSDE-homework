function [estimator,variance,conf_interval] = romberg_monte_carlo(functional, n, alpha, beta, T, dim, niveau)
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