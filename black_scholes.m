function price_series = black_scholes(dim, times, cov_matrix, interest, start_price)
    price_series = black_scholes_from_bm(times, brownian_motion(dim,times), cov_matrix, interest, start_price);
end

function price_series = black_scholes_from_bm(times, brownian_motion, cov_matrix, interest, start_price)
    % times: sorted row vector (1xT) of positive timestamps for prices
    % brownian_motion: d x T brownian motion where column i matches times(i)
    % cov_matrix: a m x m covariance matrix to use
    % interest: the interest rate
    % start_price: start_price column vector
    % :return price_series: a m x T price series is returned 
    vec = (interest-0.5*sum(cov_matrix.^2, 2));
    price_series = start_price.*exp(vec * times + cov_matrix * brownian_motion);
end
