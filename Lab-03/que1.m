function que1()
    T = 1; K = 20; r = 0.07; sigma = 0.2; delta = 0.01;

    S_max = 100;
    dx = 0.5; dt = 0.0005;

    lbc = @(tau) 0;
    rbc = @(tau) S_max - K*exp(-r*tau);

    Eps = [1e-2, 1e-4, 1e-6];

    for eps = Eps
        [U,x,tau] = FTCS(0, 0, S_max, T, dx, dt, eps, lbc, rbc, r, sigma, delta,K);
        t = flip(tau);
        
        plots(U,x,t,' FTCS', eps, dx);
    end

    for eps = Eps
        [U,x,tau] = BTCS(0, 0, S_max, T, dx, dt, eps, lbc, rbc, r, sigma, delta,K);
        t = flip(tau);
        
        plots(U,x,t,' BTCS', eps, dx);
    end

    for eps = Eps
        [U,x,tau] = Crank_Nicolson(0, 0, S_max, T, dx, dt, eps, lbc, rbc, r, sigma, delta,K);
        t = flip(tau);
        
        plots(U,x,t,' Crank-Nicolson', eps, dx);
    end
end