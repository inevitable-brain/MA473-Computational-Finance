function [U,s,t] = Crank_Nicolson(a,c,b,d,ds,dt,f1,g1,g2,r,delta,sigma_bar)
    s = a:ds:b;
    t = c:dt:d;
    
    N = length(s);
    M = length(t);
    
    U = zeros(M,N);
    U(:,1) = g1(t);
    U(:,N) = g2(t);
    U(1,:) = f1(s);

    alpha = @(s) sigma_bar(s) * s * (1 - s);
    beta = @(s) (r-delta) * s * (1 - s);
    gamma = @(s) r*(1 - s) + delta*s;
    
    c1 = @(s) (-0.5 * alpha(s)^2 * dt / ds^2 - 0.5 * beta(s) * dt / ds) / 2;
    c2 = @(s) 1 + 0.5 * (alpha(s)^2 * dt / (ds)^2 + gamma(s) * dt);
    c3 = @(s) (-0.5 * alpha(s)^2 * dt / ds^2 + 0.5 * beta(s) * dt / ds) / 2;

    d1 = @(s) (+0.5 * alpha(s)^2 * dt / ds^2 + 0.5 * beta(s) * dt / ds) / 2;
    d2 = @(s) 1 - 0.5 * (alpha(s)^2 * dt / (ds)^2 + gamma(s) * dt);
    d3 = @(s) (+0.5 * alpha(s)^2 * dt / ds^2 - 0.5 * beta(s) * dt / ds) / 2;

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
