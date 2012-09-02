function res = geman_func(R)
    res = 0.5 * (R .^ 2) ./ (1 + R .^ 2);
end