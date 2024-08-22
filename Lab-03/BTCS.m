function [U,x,t] = BTCS(a,c,b,d,dx,dt,eps,g1,g2,r,sigma,delta,K)
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x,K,eps);
    
    c1 = @(x) (-0.5*(sigma^2*x^2)/dx^2 - 0.5*(r-delta)*x/dx);
    c2 = @(x) ((1/dt)+(r)+(sigma^2*x^2)/dx^2);
    c3 = @(x) (-0.5*(sigma^2*x^2)/dx^2 + 0.5*(r-delta)*x/dx);
    
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
        
        b(2:N-1) = U(i-1,2:N-1)/dt;
        b(1) = U(i,1);
        b(end) = U(i,end);

        U(i,:) = (A\b)';
    end
end