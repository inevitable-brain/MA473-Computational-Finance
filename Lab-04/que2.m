function que2()
    T = 1; r = 0.04; delta = 0.1; q = 0.2;
    sigma_bar = @(x) q*x / (4*(1-x));

    dx = 0.005; dt = 0.005;

    % European Put
    lbc = @(tau) exp(-r*tau);
    rbc = @(tau) 0;
    f1 = @(x) max(1 - 2*x, 0);

    [U,x,tau] = FTCS(0, 0, 1, T, dx, dt, f1, lbc, rbc, r, delta,sigma_bar);
    [V,S,t] = transform(U, x, tau, q, T);
    
    plots(V,S,t,2,' FTCS for Put',dx);

    [U,x,tau] = BTCS(0, 0, 1, T, dx, dt, f1, lbc, rbc, r, delta,sigma_bar);
    [V,S,t] = transform(U, x, tau, q, T);
    
    plots(V,S,t,2,' BTCS for Put',dx);

    [U,x,tau] = Crank_Nicolson(0, 0, 1, T, dx, dt, f1, lbc, rbc, r, delta,sigma_bar);
    [V,S,t] = transform(U, x, tau, q, T);
    
    plots(V,S,t,2,' Crank-Nicolson for Put',dx);
end