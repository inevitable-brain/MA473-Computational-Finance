function [U,x,t] = Crank_Nicolson(a,c,b,d,dx,dt,alpha,f1,g1,g2)
    lambda = (alpha*dt)/(2*dx^2);
    
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t(1:M));
    U(:,N) = g2(t(1:M));
    U(1,1:N) = f1(x(1:N));

    A=zeros(N-2,N-2);
    W=zeros(N-2,1);
    A(1,1) = (1+2*lambda);
    A(1,2) = -lambda;
    A(N-2,N-3) = -lambda;
    A(N-2,N-2) = (1+2*lambda);
    
    for i=2:N-3
        A(i,i-1)=-lambda;
        A(i,i)=(1+2*lambda);
        A(i,i+1)=-lambda;
    end
    for i=1:M-1
        W(1) = (1-2*lambda)*U(i,2) + lambda*U(i,3) + lambda*(U(i,1)+U(i+1,1));
        for j=2:N-3
            W(j) = lambda*U(i,j) + (1-2*lambda)*U(i,j+1) + lambda*U(i,j+2);
        end
        W(N-2) = lambda*U(i,N-2) + (1-2*lambda)*U(i,N-1) + lambda*(U(i,N)+U(i+1,N));
        U(i+1,2:N-1) = A\W;
    end
end