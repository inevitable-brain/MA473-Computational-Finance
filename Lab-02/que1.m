function que1()
    % given values
    T = 1; K = 10; r = 0.06; sigma = 0.3; delta = 0;
    
    % taking x_min, x_max and defining tau_min, tau_max
    x_min = -2; x_max = 2;
    tau_min = 0; tau_max = sigma^2/2;
    dx = 0.05; dt = dx^2/2;
    
    q_delta = 2 * (r - delta)/sigma^2;
    
    % defining initial and boundary conditions
    init_cond = @(x) max(exp(x * (q_delta + 1) / 2) - exp(x * (q_delta - 1) / 2), 0);
    lbc = @(tau) 0;
    rbc = @(tau) exp((q_delta + 1) * x_max/2 + (q_delta + 1)^2 * tau/4);
    
    % FTCS
    [U, x, tau] = FTCS(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
    
    [v, S, t] = transform(x, tau, U, T, K, r, sigma, delta);
    
    plots(v, S, t, x, tau, ' FTCS');
    
    % BTCS
    [U, x, tau] = BTCS(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
    
    [v, S, t] = transform(x, tau, U, T, K, r, sigma, delta);
    
    plots(v, S, t, x, tau, ' BTCS');
    
    % Crank-Nicolson
    [U, x, tau] = Crank_Nicolson(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
    
    [v, S, t] = transform(x, tau, U, T, K, r, sigma, delta);
    
    plots(v, S, t, x, tau, ' Crank-Nicolson');
    
end