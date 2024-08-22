function [U,x,t] = FTCS(a,c,b,d,dx,dt,f1,g1,g2,r,delta,sigma_bar)
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x);
    
    for i=1:M-1
        for j=2:N-1
            alpha = sigma_bar(x(j)) * x(j) * (1 - x(j));
            beta = (r-delta) * x(j) * (1 - x(j));
            gamma = r*(1 - x(j)) + delta*x(j);

            a = 0.5 * alpha^2 * dt / (dx)^2 + beta * dt / (2*dx);
            b = 1 - alpha^2 * dt / (dx)^2 - gamma * dt;
            c = 0.5 * alpha^2 * dt / (dx)^2 - beta * dt / (2*dx);
            U(i+1,j) = a*U(i,j+1) + b*U(i,j) + c*U(i,j-1);
        end
    end
end