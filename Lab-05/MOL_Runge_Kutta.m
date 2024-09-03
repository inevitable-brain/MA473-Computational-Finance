function [U,x,t] = MOL_Runge_Kutta(a,c,b,d,dx,dt,alpha,f1,g1,g2)
    lambda = (alpha*dt)/(dx^2);
    
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(x);

    A = diag(-2*ones(1,N-2)) + diag(ones(1,N-3), 1) + diag(ones(1,N-3), -1);
    beta = zeros(N-2,1);
    c2 = 1; a21 = 1;
    w1 = 1/2; w2 = 1/2;
    
    for i=1:M-1
        beta(1,1) = g1(t(i));
        beta(N-2,1) = g2(t(i));
        K1 = lambda * (A*U(i,2:N-1)' + beta);

        beta(1,1) = g1(t(i) + c2*dt);
        beta(N-2,1) = g2(t(i) + c2*dt);
        K2 = lambda * (A * (U(i,2:N-1)' + a21*K1) + beta);
        
        U(i+1, 2:N-1) = U(i, 2:N-1) + w1*K1' + w2*K2';
    end
end