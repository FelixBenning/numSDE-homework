function simulation = ornsteinUhlenbeck(T, x_0, mu, lambda, n)
    times = 0:T/n:T;
    cov_matrix = covariance(times, lambda);
    A = chol(cov_matrix);
    return expectedValue(times, lambda) + rand(1, length(times)) * A;
end

function covaraince(times, lambda)
    % times: row vector of times (1xT)
    % return: TxT covariance matrix corresponding to times
        
end

function expectedValue(times, lambda)
    % times: row vector of times
    % return: row vector of the expected value corresponding to times

end