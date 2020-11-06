function simulation = ornsteinUhlenbeck(T, x_0, mu, lambda, n)
    times = 0:T/n:T;
    cov_matrix = covariance(times, lambda);
    A = chol(cov_matrix);
    simulation = expectedValue(times, lambda) + rand(1, length(times)) * A;
end

function cov_matrix = covaraince(times, lambda)
    % times: row vector of times (1xT)
    % return: TxT covariance matrix corresponding to times
    
    exp_vec = exp(-lambda .* times);
    cov_matrix = transpose(exp_vec) * exp_vec ./(2*lambda);
    
    temp = repelem( transpose(exp(2*lambda .* times)),1, length(times));
    min_matrix = triu(temp) + transpose(triu(temp,1)) - 1;
    
    cov_matrix = cov_matrix .* min_matrix;
end

function expectation = expectedValue(times, lambda)
    % times: row vector of times
    % return: row vector of the expected value corresponding to times
    
end