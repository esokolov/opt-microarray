function L = huber_func(x)
    L = x;
    idx = abs(L) <= 1;
    L(idx) = 0.5 * (L(idx) .^ 2);
    L(~idx) = abs(L(~idx)) - 0.5;
end