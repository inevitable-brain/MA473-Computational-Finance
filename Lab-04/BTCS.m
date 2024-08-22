function [U,x,t] = BTCS(a,c,b,d,dx,dt,f1,g1,g2,r,delta,sigma_bar)
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x);

    alpha = @(x) sigma_bar(x) * x * (1 - x);
    beta = @(x) (r-delta) * x * (1 - x);
    gamma = @(x) r*(1 - x) + delta*x;
    
    c1 = @(x) -0.5 * alpha(x)^2 * dt / dx^2 - 0.5 * beta(x) * dt / dx;
    c2 = @(x) 1 + alpha(x)^2 * dt / (dx)^2 + gamma(x) * dt;
    c3 = @(x) -0.5 * alpha(x)^2 * dt / dx^2 + 0.5 * beta(x) * dt / dx;
    
    for i = 2:M
        A = zeros(N, N);
        b = zeros(N, 1);
        
        for j = 2:N-1
            A(j,j+1) = c1(x(j));
            A(j,j) = c2(x(j));
            A(j,j-1) = c3(x(j));
        end
        A(1,1) = 1; A(1,2) = 0;
        A(N,N) = 1; A(N,N-1) = 0;
        
        b(2:N-1) = U(i-1,2:N-1);
        b(1) = U(i,1);
        b(end) = U(i,end);

        U(i,:) = (A\b)';
    end
end