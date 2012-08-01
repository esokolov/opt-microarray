z = 75;
arr_idx=484;
alpha = -0.5;
beta  = -0.5;
eps  =  1e-12;
F = @(z) (A * z) ./ (1 + B * z);
f = @(z) sum(sum(- (1 / (alpha * beta)) * ...
            (bsxfun(@times, bsxfun(@power, I(:,arr_idx) + eps, alpha), bsxfun(@power, F(z) + eps, beta)) - ...
            (alpha / (alpha + beta)) * bsxfun(@power, I(:,arr_idx) + eps, alpha + beta) - ...
            (beta / (alpha + beta)) * bsxfun(@power, F(z) + eps, alpha + beta))));
fd = @(z) (1 ./ (alpha * z .^ 2)) .* sum(bsxfun(@times, 1 ./ A, (F(z) .^ (beta + 1)) .* ...
                (F(z) .^ alpha - (I(:,arr_idx) + eps) .^ alpha)), 1);
sd = @(z) (1 ./ (alpha * z .^ 2)) .* sum((1 ./ ((1 + B * z) .^ 2)) .* (F(z) .^ beta) .* ...
                ((F(z) .^ alpha) .* (alpha - 2 * B * z + beta - 1) + ...
                ((I(:,arr_idx) + eps) .^ alpha) .* (2 * B * z - beta + 1)), 1);
a = @(z) -0.5 * sd(z) / fd(z);
b = @(z) - fd(z)^2 / (sd(z) * exp(sd(z)/fd(z)*z^2));
d = @(z) f(z) - fd(z)^2/sd(z);
approx = @(c,z) d(z) - b(z) * exp(-a(z)*c^2);

c=10:1:400;
y1 = zeros(size(c));
y2 = zeros(size(c));
for i=1:length(c)
    y1(i) = f(c(i));
    y2(i) = approx(c(i),z);
end
figure
plot(c,y1);
hold on
plot(c,y2,'r');
