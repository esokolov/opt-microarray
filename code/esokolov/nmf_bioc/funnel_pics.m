array_idx = 5;
func = @(c) nmf_alpha_beta_divergence(I(:, array_idx), A * c ./ (1 + B * c), alpha, beta);
func_grad = @(c) (1 ./ (alpha * c .^ 2)) .* sum(bsxfun(@times, 1 ./ A, ((A * c ./ (1 + B * c)) .^ (beta + 1)) .* ...
    ((A * c ./ (1 + B * c)) .^ alpha - (I(:, array_idx) + eps) .^ alpha)), 1);
func_hess = @(c) (1 ./ (alpha * c .^ 2)) .* sum((1 ./ ((1 + B * c) .^ 2)) .* ...
    ((A * c ./ (1 + B * c)) .^ beta) .* (((A * c ./ (1 + B * c)) .^ alpha) .* (alpha - 2 * B * c + beta - 1) + ...
                ((I(:, array_idx) + eps) .^ alpha) .* (2 * B * c - beta + 1)), 1);
            
c_curr = C(array_idx);
c_curr = 500;

%a = - (- func_grad(c_curr)^2 / func(c_curr) + func_hess(c_curr)) / (2 * func(c_curr));
%b = -2 * a * c_curr - func_grad(c_curr) / func(c_curr);
%c_coef = -log(-func(c_curr) * exp(a * c_curr^2 + b * c_curr));
%c_coef = 0;

%a = 0.5 * ((func_grad(c_curr) / func(c_curr)) ^ 2 - func_hess(c_curr) / func(c_curr));
%b = -2 * a * c_curr - func_grad(c_curr) / func(c_curr);
%c_coeg = -log(-

func = @(x) x^2;
func_grad = @(x) 2*x;
func_hess = @(x) 2;

z1 = -2;
z2 = 2;

g = func(z1) - func(z2);
gamma1 = roots([(func_hess(z2) - func_hess(z1)); (func_hess(z2) * g + func_grad(z2)^2 - func_grad(z1)^2 - 2 * g * func_hess(z1)); ...
    (-(g^2 * func_hess(z1) - 2 * g * func_grad(z1)^2)); (-g^2 * func_grad(z1)^2)]);
gamma1 = 10;
gamma2 = func(z1) - func(z2) + gamma1;
idx = find(gamma1 > 0 & gamma2 > 0, 1, 'first');
%gamma = max(gamma1(idx), gamma2(idx));
if (gamma1(idx) > gamma2(idx))
    fval = func(z1);
    gval = func_grad(z1);
    hval = func_hess(z1);
    gamma = gamma1(idx);
    z = z1;
else
    fval = func(z2);
    gval = func_grad(z2);
    hval = func_hess(z2);
    gamma = gamma2(idx);
    z = z2;
end

a = (gamma * hval + gval^2) / (gamma^2);
b = gval / gamma - a * z;
new_z = - b / a;
f0 = fval + gamma;
c_val = 0.5 * a * new_z ^2 + b * new_z;
gamma_coef = (f0 - fval) * exp(0.5 * a * z ^2 + b * z + c_val);


func_approx = @(c) f0 - gamma_coef * exp(-0.5 * a * c^2 - b * c - c_val);

x = -5:0.01:5;
y = x;
y1 = x;
for j = 1:length(x)
    y(j) = func(x(j));
    y1(j) = func_approx(x(j));
end
plot(x, y);
hold on;
plot(x, y1, 'r');
%plot([c_curr, c_curr], [0 0.5], '--');