function [c0, c1, c2, c4, c6, c8] = pi_epsilon(eps)
    c0 = 35 * eps/256;
    c1 = 1 / 2;
    c2 = 35 / (64*eps);
    c4 = -35 / (128*eps^3);
    c6 = 7 / (64*eps^5);
    c8 = -5 / (256*eps^7);
end