function simulation = ornsteinUhlenbeck(T, x_0, mu, lambda, n)
    times = T/n:T/n:T;
    [expected_value, cov_matrix] = distribution_parameters(times, lambda, mu, x_0);
    A = chol(cov_matrix);
    simulation = [x_0, expected_value + randn(1, length(times)) * A];
end

function [expected_value, cov_matrix] = distribution_parameters(times, lambda, mu, x_0)
    % times: row vector of times (1xT)
    % return: TxT covariance matrix corresponding to times
    
    exp_vec = exp(-lambda .* times);
    
    % cov matrix calculation
    cov_matrix = transpose(exp_vec) * exp_vec ./(2*lambda);
    
    temp = repelem( transpose(exp(2*lambda .* times)), 1, length(times));
    min_matrix = triu(temp) + transpose(triu(temp,1)) - 1;
    
    cov_matrix = cov_matrix .* min_matrix;
    
    % expected_value calculation
    
    expected_value = x_0 .* exp_vec + (mu .*(1 - exp_vec));
    
end