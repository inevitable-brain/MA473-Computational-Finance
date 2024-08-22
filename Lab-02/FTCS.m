function [U,x,t] = FTCS(a,c,b,d,dx,dt,alpha,f1,g1,g2)
    lambda = (alpha*dt)/(dx^2);
    
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
            U(i+1,j) = lambda*U(i,j+1) + (1-2*lambda)*U(i,j) + lambda*U(i,j-1);
        end
    end
end