function que1()
    % given values
    T = 1; K = 10; r = 0.06; sigma = 0.3; delta = 0;
    
    % taking x_min, x_max and defining tau_min, tau_max
    x_min = -3; x_max = 3;
    tau_min = 0; tau_max = sigma^2/2;

    dx_array = 0.03:0.002:0.1; dt_array = dx_array / 140;

    dx = 0.03; dt = dx/200;
    
    q_delta = 2 * (r - delta)/sigma^2;
    
    % defining initial and boundary conditions
    init_cond = @(x) max(exp(x * (q_delta + 1) / 2) - exp(x * (q_delta - 1) / 2), 0);
    lbc = @(tau) 0;
    rbc = @(tau) exp((q_delta + 1) * x_max/2 + (q_delta + 1)^2 * tau/4);
    
    % MOL Explicit Euler
    [U_ex, x, tau] = MOL_Explicit_Euler(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);

    [v_ex, S, t] = transform(x, tau, U_ex, T, K, r, sigma, delta);

    % plots(v_ex, S, t, ' MOL Explicit-Euler');

    % MOL Implicit Euler
    [U_im, x, tau] = MOL_Implicit_Euler(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
    
    [v_im, S, t] = transform(x, tau, U_im, T, K, r, sigma, delta);
    
    % plots(v_im, S, t, ' MOL Implicit-Euler');

    % MOL Runge-Kutta
    [U_RK, x, tau] = MOL_Runge_Kutta(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
    
    [v_RK, S, t] = transform(x, tau, U_RK, T, K, r, sigma, delta);
    
    % plots(v_RK, S, t, ' MOL Runge-Kutta');

    % Plotting Order Of Convergence
    E = zeros(1, length(dx_array));
    E2 = zeros(1, length(dx_array));

    fun = @(x,tau) (exp((q_delta - 1) * (2*x + (q_delta - 1) * tau)) * (exp(x + q_delta*tau)*(erf((x + (q_delta + 1) * tau) / (2*sqrt(tau))) + 1) - erf((x + (q_delta - 1) * tau) / (2*sqrt(tau))) - 1 )) / 2;
    for i = 1:length(dx_array)
        dx = dx_array(i);
        dt = dt_array(i);

        dx2 = dx_array(i) / 2;

        [U, x_, tau] = MOL_Explicit_Euler(x_min, tau_min, x_max, tau_max, dx, dt, 1, init_cond, lbc, rbc);
        [U2, x2, tau2] = MOL_Explicit_Euler(x_min, tau_min, x_max, tau_max, dx2, dt, 1, init_cond, lbc, rbc);

        U_actual = zeros(length(tau), length(x_));
        for j = 1:length(x_)
            for k = 1:length(tau)
                U_actual(k, j) = fun(x_(j), tau(k));
            end
        end

        U_actual2 = zeros(length(tau2), length(x2));
        for j = 1:length(x2)
            for k = 1:length(tau2)
                U_actual2(k, j) = fun(x2(j), tau2(k));
            end
        end

        % [v, ~, ~] = transform(x, tau, U, T, K, r, sigma, delta);
        % [v_actual, ~, ~] = transform(x, tau, U_actual, T, K, r, sigma, delta);
        % [x_mesh, tau_mesh] = meshgrid(x_,tau);
        % 
        % figure;
        % mesh(x_mesh, tau_mesh, U_actual,'FaceColor','flat', 'FaceAlpha',0.6);
        % hold off;

        E(1,i) = max(max(abs(U - U_actual)));
        E2(1,i) = max(max(abs(U2 - U_actual2)));
    end

    figure;
    plot(dx_array, log2(E./E2), 'LineWidth',1);
    hold off;
end