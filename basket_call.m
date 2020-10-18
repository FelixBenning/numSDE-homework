
function payoff = basket_call(price_at_expiration, weights, strike, expiration, interest)
    % price_at_expiration: (dx1) price vector of the d products in the
    % basket
    % weights: (1xd) weights of the d products in the basket
    % strike: cost of the entire basket fixed at time 0
    % expiration: expiration time
    % interest: interest rate
    % :return payoff: payoff of the basket call
    payoff = exp(-interest*expiration)*(weights*price_at_expiration - strike);
    payoff = max(payoff, 0);
end