 function [U,x,t] = MOL_Explicit_Euler(a,c,b,d,dx,dt,alpha,f1,g1,g2)
    lambda = (alpha*dt)/(dx^2);
    
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x);

    I = diag(ones(1,N-2));
    A = diag(-2*ones(1,N-2)) + diag(ones(1,N-3), 1) + diag(ones(1,N-3), -1);
    beta = zeros(N-2,1);
    
    for i=1:M-1
        beta(1,1) = g1(t(i));
        beta(N-2,1) = g2(t(i));
        U(i+1, 2:N-1) = ((I + lambda*A) * U(i,2:N-1)' + lambda*beta)';
    end
end