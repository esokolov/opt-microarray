if (alpha == 0 && beta == 0)
        
elseif (alpha == 0)

elseif (beta == 0)

elseif (alpha == -beta)

else
    syms t a b c f g h alpha_sym beta_sym inten_sym;
    langmuir_func_sym = symfun(sym('(a + t * f) * (c + t * h) / (1 + (b + t * g) * (c + t * h))'), [t a b c f g h]);
    %divergence_summand = symfun(sym(strcat('-(1 / (alpha_sym * beta_sym)) * ', ...
    %    '(inten_sym^alpha * langmuir_func_sym(t, a, b, c, f, g, h)^beta - (alpha_sym / (alpha_sym + beta_sym)) * inten_sym ^ (alpha_sym + beta_sym) - ', ...
    %    '(beta_sym / (alpha_sym + beta_sym)) * langmuir_func_sym(t, a, b, c, f, g, h) ^ (alpha_sym + beta_sym))')), [t a b c f g h inten_sym alpha_sym beta_sym]);
    divergence_summand = symfun(sym(strcat('-(1 / (alpha_sym * beta_sym)) * ', ...
       '(inten_sym^alpha_sym * ((a + t * f) * (c + t * h) / (1 + (b + t * g) * (c + t * h)))^beta_sym - (alpha_sym / (alpha_sym + beta_sym)) * inten_sym ^ (alpha_sym + beta_sym) - ', ...
       '(beta_sym / (alpha_sym + beta_sym)) * ((a + t * f) * (c + t * h) / (1 + (b + t * g) * (c + t * h))) ^ (alpha_sym + beta_sym))')), [t a b c f g h inten_sym alpha_sym beta_sym]);
    divergence_diff_summand = diff(divergence_summand, t);
    matlabFunction(divergence_diff_summand, 'file', 'divergence_diff_summand_code.m', 'vars', {'t', 'a', 'b', 'c', 'f', 'g', 'h', 'inten_sym', 'alpha_sym', 'beta_sym'});
    matlabFunction(divergence_summand, 'file', 'divergence_summand_code.m', 'vars', {'t', 'a', 'b', 'c', 'f', 'g', 'h', 'inten_sym', 'alpha_sym', 'beta_sym'});
end