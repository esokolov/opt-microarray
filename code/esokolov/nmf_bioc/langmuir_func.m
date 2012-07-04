function res = langmuir_func(A, B, C)
    res = (A * C) ./ (1 + B * C);
end