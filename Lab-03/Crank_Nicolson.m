function [U,s,t] = Crank_Nicolson(a,c,b,d,ds,dt,eps,g1,g2,r,sigma,delta,K)
    s = a:ds:b;
    t = c:dt:d;

    N = length(s);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(s,K,eps);

    c1 = @(s) (-0.5*(sigma^2*s^2)/ds^2 - 0.5*(r-delta)*s/ds)/2;
    c2 = @(s) ((1/dt)+ 0.5*((r)+(sigma^2*s^2)/ds^2));
    c3 = @(s) (-0.5*(sigma^2*s^2)/ds^2 + 0.5*(r-delta)*s/ds)/2;

    d1 = @(s) (+0.5*(sigma^2*s^2)/ds^2 + 0.5*(r-delta)*s/ds)/2;
    d2 = @(s) ((1/dt)- 0.5*((r)+(sigma^2*s^2)/ds^2));
    d3 = @(s) (+0.5*(sigma^2*s^2)/ds^2 - 0.5*(r-delta)*s/ds)/2;

    for i=2:M
        A=zeros(N,N);
        W=zeros(N,1);
        A(1,1) = 1; A(1,2) = 0;
        A(N,N) = 1; A(N,N-1) = 0;
        
        for j=2:N-1
            A(j,j+1) = c1(s(j));
            A(j,j) = c2(s(j));
            A(j,j-1) = c3(s(j));
        end

        W(1) = U(i,1);
        W(end) = U(i,end);
        for j=2:N-1
            W(j) = d2(s(j))*U(i-1,j) + d1(s(j))*U(i-1,j+1) + d3(s(j))*U(i-1,j-1);
        end
        U(i,:) = (A\W)';
    end
end
