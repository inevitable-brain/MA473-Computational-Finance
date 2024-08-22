function que1()
    fun = @(x) cos(pi*x/2); g1 = @(t) t-t; g2 = @(t) t-t;

    actual_f = @(x,t) exp(-1*(pi^2).*t).*cos(pi*x/2);
    h_vals = [1e-2, 1e-3, 1e-4];
    k_vals = [5e-4, 1e-3, 1e-2];
    
    % for FTCS matlab limits are huge with the three combinations so we
    % check with other values
    % 
    [U,x,t] = FTCS(-1,0,2,1,0.2,5e-3,4,fun,g1,g2);
    plots(x,t,U,actual_f,'Numerical and exact plot using FTCS');

    [x,t,U] = FTCS(-1,0,2,1,0.1,0.1,4,fun,g1,g2);
    plots(x,t,U,actual_f,'Numerical and exact plot using FTCS');

    [x,t,U] = FTCS(-1,0,2,1,0.1,0.5,4,fun,g1,g2);
    plots(x,t,U,actual_f,'Numerical and exact plot using FTCS');

    for i = 1:3
        h = h_vals(i);
        k = k_vals(i);

        [U,x,t] = BTCS(-1,0,2,1,h,k,4,fun,g1,g2);
        plots(x,t,U,actual_f,'Numerical and exact plot using BTCS');

        [U,x,t] = Crank_Nicolson(-1,0,2,1,h,k,4,fun,g1,g2);
        plots(x,t,U,actual_f,'Numerical and exact plot using Crank Nicolson');

    end
    % 
end



