function terminal_cond = f1(x, K, eps)
    terminal_cond = zeros(1,length(x));

    [c0, c1, c2, c4, c6, c8] = pi_epsilon(eps);

    for i = 1:length(x)
        if x(i)-K >= eps
            terminal_cond(i) = x(i) - K;
        elseif x(i)-K <= eps
            terminal_cond(i) = 0;
        else
            terminal_cond(i) = c0 + c1*(x(i)-K) + c2*(x(i)-K).^2 + c4*(x(i)-K).^4 + c6*(x(i)-K).^6 + c8*(x(i)-K).^8;
        end
    end
end