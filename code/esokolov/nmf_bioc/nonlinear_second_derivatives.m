%%
syms a b c alpha_sym beta_sym inten_sym;
%divergence_summand = symfun(sym(strcat('-(1 / (alpha_sym * beta_sym)) * ', ...
%    '(inten_sym^alpha * langmuir_func_sym(t, a, b, c, f, g, h)^beta - (alpha_sym / (alpha_sym + beta_sym)) * inten_sym ^ (alpha_sym + beta_sym) - ', ...
%    '(beta_sym / (alpha_sym + beta_sym)) * langmuir_func_sym(t, a, b, c, f, g, h) ^ (alpha_sym + beta_sym))')), [t a b c f g h inten_sym alpha_sym beta_sym]);
divergence_summand = symfun(sym(strcat('-(1 / (alpha_sym * beta_sym)) * ', ...
   '(inten_sym^alpha_sym * (a * c / (1 + b * c))^beta_sym - (alpha_sym / (alpha_sym + beta_sym)) * inten_sym ^ (alpha_sym + beta_sym) - ', ...
   '(beta_sym / (alpha_sym + beta_sym)) * (a * c / (1 + b * c)) ^ (alpha_sym + beta_sym)) + 1e-5 * (1 + b * c)^2 / a')), [a b c inten_sym alpha_sym beta_sym]);
%divergence_diff_summand = diff(divergence_summand, t);
%matlabFunction(divergence_diff_summand, 'file', 'divergence_diff_summand_code.m', 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
%matlabFunction(divergence_summand, 'file', 'divergence_summand_code.m', 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});

da = matlabFunction(diff(divergence_summand, a), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
db = matlabFunction(diff(divergence_summand, b), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dc = matlabFunction(diff(divergence_summand, c), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});

daa = matlabFunction(diff(diff(divergence_summand, a), a), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dab = matlabFunction(diff(diff(divergence_summand, a), b), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dac = matlabFunction(diff(diff(divergence_summand, a), c), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dba = matlabFunction(diff(diff(divergence_summand, b), a), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dbb = matlabFunction(diff(diff(divergence_summand, b), b), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dbc = matlabFunction(diff(diff(divergence_summand, b), c), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dca = matlabFunction(diff(diff(divergence_summand, c), a), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dcb = matlabFunction(diff(diff(divergence_summand, c), b), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});
dcc = matlabFunction(diff(diff(divergence_summand, c), c), 'vars', {'a', 'b', 'c', 'inten_sym', 'alpha_sym', 'beta_sym'});

%%
alpha = -0.5;
beta = -0.5;

Abig = repmat(A, [1, size(I, 2)]);
Bbig = repmat(B, [1, size(I, 2)]);
Cbig = repmat(C, [size(I, 1), 1]);

grad = [sum(sum(da(Abig, Bbig, Cbig, I, alpha, beta)));
        sum(sum(db(Abig, Bbig, Cbig, I, alpha, beta)));
        sum(sum(dc(Abig, Bbig, Cbig, I, alpha, beta))) ];

hess = [sum(sum(daa(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dab(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dac(Abig, Bbig, Cbig, I, alpha, beta)));
        sum(sum(dba(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dbb(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dbc(Abig, Bbig, Cbig, I, alpha, beta)));
        sum(sum(dca(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dcb(Abig, Bbig, Cbig, I, alpha, beta))), sum(sum(dcc(Abig, Bbig, Cbig, I, alpha, beta))) ];
    
%%
grad_a = size(A);
grad_b = size(B);
grad_c = size(C);
for i = 1:length(A)
    grad_a(i) = sum(da(A(i), B(i), C, I(i, :), alpha, beta));
    grad_b(i) = sum(db(A(i), B(i), C, I(i, :), alpha, beta));
end
for i = 1:length(C)
    grad_c(i) = sum(dc(A, B, C(i), I(:, i), alpha, beta));
end