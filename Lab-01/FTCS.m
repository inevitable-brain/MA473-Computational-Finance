function [U,x,t] = FTCS(a,c,b,d,dx,dt,alpha,f1,g1,g2)
    N = round(b/dx + 1);
    M = round(d/dt + 1);
    lambda = (alpha*dt)/(dx^2);
    
    x = zeros(N,1);
    t = zeros(M,1);
    for i=1:N
        x(i) = a + (i-1)*dx;
    end
    for i=1:M
        t(i) = c + (i-1)*dt;
    end
    
    U = zeros(M,N);
    U(:,1) = g1(t(1:M));
    U(:,N) = g2(t(1:M));
    U(1,1:N) = f1(x(1:N));
    
    for i=1:M-1
        for j=2:N-1
            U(i+1,j) = lambda*U(i,j+1) + (1-2*lambda)*U(i,j) + lambda*U(i,j-1);
        end
    end
end