function [x,v] = leapfrog(fun, delta, x_0, v_0, T, eta, sigma)

    update_fun = @(nu, x, v, dBM) v+0.5*((fun(x+0.5*nu)-eta*nu)*delta+sigma*dBM);
    implicit_fun = @(nu, x, v, dBM) update_fun(nu,x,v,dBM) - nu;
    
     
    times = 0:delta:T;
    bm = brownian_motion(1,times);
    d_bm = diff(bm);
    
    x = zeros(1, length(times));
    v = zeros(1, length(times));
    x(1) = x_0;
    v(1) = v_0;
    for idx = 1:(length(times)-1)
        nu = fzero(@(nu) implicit_fun(nu, x(idx), v(idx), d_bm(idx)), 0);
        x(idx+1) = x(idx)+ nu*delta;
        v(idx+1) = update_fun(nu, x(idx), v(idx), d_bm(idx));
    end
    
    
end