function [V, S, t] = transform(V_bar, sai, tau, q, T)
    S = q*sai ./ (1 - sai);
    t = T - tau;
    V = (S + q) .* V_bar(:, :);
end