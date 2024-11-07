function [U, x, t] = MOHL_Runge_Kutta1(a,c,b,d,dx,dt,f1,g1,g2,q_delta)
    x = a:dx:b;
    t = c:dt:d;
    
    N = length(x);
    M = length(t);

    U = zeros(M, N);
    maxIter=100;
    U(1:end, 1) = g1(t);
    U(1:end, end) = g2(t);
    U(1, 1:end) = f1(x);

    z = 0;
    for j=2:N-1
        if x(j) <= 0
            z = 0;
        else
            z = 0.5*(q_delta+1)*exp(x(j)*0.5*(q_delta+1)) - 0.5*(q_delta-1)*exp(x(j)*0.5*(q_delta-1));
        end
        for i = 2:M
            K1 = dx*z;
            l1 = (dx/dt)*(U(i, j-1) - U(i-1, j-1));
            K2 = dx*(z + l1);
            l2 = (dx/dt)*(U(i, j-1) - U(i-1, j-1) + K1);
            U(i,j) = U(i-1, j) + K1/2 + K2/2;
            for k=1:maxIter
                z1 = z + l1/2 + l2/2;
                K1 = dx*z1;
                l1 = (dx/dt)*(U(i, j) - U(i-1, j));
                K2 = dx*(z1 + l1);
                l2 = (dx/dt)*(U(i, j) - U(i-1, j) + K1);
                U(i,j) = U(i-1, j) + K1/2 + K2/2;
            end
            z = z + l1/2 + l2/2;
        end
    end
end
