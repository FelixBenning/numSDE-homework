function payoff = G(price_series, expiration, interest)
    % price_series: (2xT) dimensional price series
    % expiration: expiration time
    % interest: interest rate
    % :return payoff: payoff of this payoff function
    payoff = exp(-interest*expiration)*(max(price_series(1,:))-min(price_series(2,:)));
    payoff = max(payoff, 0);
end