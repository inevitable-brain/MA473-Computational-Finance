function [U,x,t] = FTCS(a,c,b,d,dx,dt,eps,g1,g2,r,sigma,delta,K)
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x,K,eps);
    
    for i=1:M-1
        for j=2:N-1
            a = 0.5 * sigma^2 * x(j)^2 * dt / (dx)^2 + (r-delta) * x(j) * dt / (2*dx);
            b = 1 - sigma^2 * x(j)^2 * dt / (dx)^2 - r * dt;
            c = 0.5 * sigma^2 * x(j)^2 * dt / (dx)^2 - (r-delta) * x(j) * dt / (2*dx);
            U(i+1,j) = a*U(i,j+1) + b*U(i,j) + c*U(i,j-1);
        end
    end
end