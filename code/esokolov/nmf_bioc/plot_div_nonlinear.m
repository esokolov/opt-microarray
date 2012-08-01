%% plot divergence(B)
%f1 = @(b) sum(log(I(16, :)) -1 - log(A(16) * C ./ (1 + b * C)) + A(16) * C ./ (I(16, :) .* (1 + b * C)));
array_idx = 4;
f_b_plot = @(b) nmf_alpha_beta_divergence(I(array_idx, :), A(array_idx) * C ./ (1 + b * C), alpha, beta);
x = 1:0.01:100;
y = x;
for j = 1:length(x)
    y(j) = f_b_plot(x(j));
end
figure; plot(x, y)

%% plot divergence(C)
probe_idx = 1;
f_c_plot = @(c) nmf_alpha_beta_divergence(I(probe_idx, :), A * c ./ (1 + B * c), alpha, beta);
x = 0:1:1000;
y = x;
for j = 1:length(x)
    y(j) = f_c_plot(x(j));
end
figure; plot(x, y)

%% C quad approx
array_idx = 720;
%array_idx = 57;
%array_idx = 80;
%array_idx = 700;
%array_idx = 584;
f_c_plot = @(c) sum(A * c ./ (I(:, array_idx) .* (1 + B * c)) - log(A * c ./ (1 + B * c)));
f_c_reg_plot = @(c) sum(A * c ./ (I(:, array_idx) .* (1 + B * c)) - log(A * c ./ (1 + B * c))) + 1e-9 * c ^ 2;
f_c_approx = @(c, c0) sum(A * c0 ./ (I(:, array_idx) .* (1 + B * c0)) + (A ./ (I(:, array_idx) .* (1 + B * c0) .^ 2)) * (c - c0) - ...
    (A .* B ./ (I(:, array_idx) .* (1 + B * c0) .^ 3)) * ((c - c0) .^ 2) - ...
    log(A * c0 ./ (1 + B * c0)) - (1 ./ (c0 + B * (c0^2))) * (c - c0) + ...
    ((B * c0 + 0.5) ./ (c0^2 * (1 + B * c0) .^ 2)) * ((c - c0) .^ 2));
x = 1e10:1e17:1e19;
y = x;
y1 = x;
y_reg = x;
for j = 1:length(x)
    y(j) = f_c_plot(x(j));
    y_reg(j) = f_c_reg_plot(x(j));
    y1(j) = f_c_approx(x(j), C(array_idx));
end
figure;
plot(x, y, 'r', 'LineWidth', 2);
hold on;
%plot(x, y_reg, 'b', 'LineWidth', 2);
plot(x, y1, 'LineWidth', 2);
%plot([C(array_idx) C(array_idx)], [min(min(y)) max(max(y))], '--', 'LineWidth', 2);
grid on;
%legend('Functional', 'Quadratic approximation');

%% C quad approx (general case)
array_idx = 912;
%array_idx = 57;
%array_idx = 80;
%array_idx = 700;
%array_idx = 584;
f_c_plot = @(c) nmf_alpha_beta_divergence(I(:, array_idx), A * c ./ (1 + B * c), alpha, beta);
x = 1000:10:20000;
y = x;
for j = 1:length(x)
    y(j) = f_c_plot(x(j));
end
figure;
plot(x, y, 'r', 'LineWidth', 2);
grid on;

%%
array_idx = 1;
f_b_conv_plot = @(b) sum(A(array_idx) * C ./ (I(array_idx, :) .* (1 + b * C)));
f_b_conc_plot = @(b) sum(log(1 + b * C));
f_b_full_plot = @(b) f_b_conv_plot(b) + f_b_conc_plot(b);
x = 0:0.000001:0.001;
y = x;
y1 = x;
y2 = x;
for j = 1:length(x)
    y(j) = f_b_full_plot(x(j));
    y1(j) = f_b_conv_plot(x(j));
    y2(j) = f_b_conc_plot(x(j));
end
figure;
plot(x, y, 'r');
hold on;
plot(x, y1);
plot(x, y2);

%%
array_idx = 1;
f_b_conv_plot = @(b) sum(A(array_idx) * C ./ (I(array_idx, :) .* (1 + b * C)));
f_b_conc_plot = @(b) sum(log(1 + b * C));
f_b_full_plot = @(b) f_b_conv_plot(b) + f_b_conc_plot(b);
%f_b_approx = @(b) sum(b .^ 2 - 2 * (((A(array_idx) * C(1) ./ (I(array_idx, 1))) - 1) ./ C(1)) * b + (A(array_idx) * C(1) ./ (I(array_idx, 1))));
rpk = (A(array_idx) * C ./ (I(array_idx, :)));
b_coef = (rpk - log(rpk) - 1) .* C ./ (1 - rpk);
f_b_approx = @(b) sum((b_coef .* C ./ (1 - rpk)) * b .^ 2 + 2 * b_coef * b + rpk);
%f_b_approx = @(b) b .^2 - 2 * (((A(array_idx) * C(1) ./ (I(array_idx, 1))) - 1) ./ C(1)) * b + 1 + log(A(array_idx) * C(1) ./ (I(array_idx, 1))) + (((A(array_idx) * C(1) ./ (I(array_idx, 1))) - 1) ./ C(1)) .^ 2;
x = 0:0.00001:0.001;
y = x;
y1 = x;
for j = 1:length(x)
    y(j) = f_b_full_plot(x(j));
    y1(j) = f_b_approx(x(j));
end
figure;
plot(x, y, 'r');
hold on;
plot(x, y1);

%%
array_idx = 1;

f_b_conv_plot = @(b) sum(A(array_idx) * C ./ (I(array_idx, :) .* (1 + b * C)));
f_b_conc_plot = @(b) sum(log(1 + b * C));
f_b_full_plot = @(b) f_b_conv_plot(b) + f_b_conc_plot(b);

rpk = (A(array_idx) * C ./ (I(array_idx, :)));
%f_b_approx = @(b) sum(rpk ./ (1 + B_old(array_idx) * C) - (rpk .* C ./ (1 + B_old(array_idx) * C) .^ 2) * (b - B_old(array_idx)) + ...
%    ((rpk .* C .^ 2) ./ (1 + B_old(array_idx) * C) .^ 3) * (b - B_old(array_idx)) .^ 2 + ...
%    log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)) ...
%    - 0.5 * (C .^ 2 ./ (1 + C * B_old(array_idx)) .^ 2) * (b - B_old(array_idx)) .^ 2);
B_old = B;
Q = langmuir_func(A, B_old, C);
f_b_approx = @(b) nmf_alpha_beta_divergence(I(array_idx, :), Q(array_idx, :), alpha, beta) + ...
    (1 / alpha) * sum(((eps + bsxfun(@rdivide, I(array_idx, :), A(array_idx))) .^ alpha) .* ((eps + bsxfun(@rdivide, 1 + B_old(array_idx) * C, C)) .^ (alpha - 1)) - ...
                bsxfun(@rdivide, C, 1 + B_old(array_idx) * C), 2) * (b - B_old(array_idx)) + ...
                (1 / alpha) * sum((alpha - 1) * (((eps + bsxfun(@rdivide, I(array_idx, :), A(array_idx))) .^ alpha) .* ((eps + bsxfun(@rdivide, 1 + B_old(array_idx) * C, C)) .^ (alpha - 2))) + ...
                bsxfun(@rdivide, C, 1 + B_old(array_idx) * C) .^ 2, 2) * ((b - B_old(array_idx)) .^ 2);

f_b_lin = @(b) sum(rpk ./ (1 + b * C) + log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)));

x = 1:1:10;
y = x;
y1 = x;
y2 = x;
for j = 1:length(x)
    y(j) = f_b_full_plot(x(j));
    y1(j) = f_b_approx(x(j));
    y2(j) = f_b_lin(x(j));
end
figure;
plot(x, y, 'r');
hold on;
plot(x, y1);
%plot(x, y2, 'b');

%% A quad approx
probe_idx = 18;

Q = langmuir_func(A, B, C);
f_a_plot = @(a) nmf_alpha_beta_divergence(I(probe_idx, :), (a * C) ./ (1 + B(probe_idx) * C), alpha, beta) + 0.5*alpha_A/2 * a^2;

%rpk = (A(array_idx) * C ./ (I(array_idx, :)));
%f_b_approx = @(b) sum(rpk ./ (1 + B_old(array_idx) * C) - (rpk .* C ./ (1 + B_old(array_idx) * C) .^ 2) * (b - B_old(array_idx)) + ...
%    ((rpk .* C .^ 2) ./ (1 + B_old(array_idx) * C) .^ 3) * (b - B_old(array_idx)) .^ 2 + ...
%    log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)) ...
%    - 0.5 * (C .^ 2 ./ (1 + C * B_old(array_idx)) .^ 2) * (b - B_old(array_idx)) .^ 2);
F1 = (A(probe_idx) * C) ./ (1 + B(probe_idx) * C);
f_a_approx = @(a) nmf_alpha_beta_divergence(I(probe_idx, :), Q(probe_idx, :), alpha, beta) + 0.5*alpha_A/2 * a^2 + ...
    (alpha_A * A(probe_idx) + (1/alpha) * (1 ./ (A(probe_idx) + eps)) .* sum(((F1 + eps) .^ beta) .* ((F1 + eps) .^ alpha - (I(probe_idx, :) + eps) .^ alpha), 2)) * (a - A(probe_idx)) + ...
                0.5 * (alpha_A + (1 / alpha) * ((1 ./ (A(probe_idx) + eps)) .^ 2) .* ...
                sum(((F1 + eps) .^ beta) .* ((alpha + beta - 1) * ((F1 + eps) .^ alpha) - (beta - 1) * ((I(probe_idx, :) + eps) .^ alpha)), 2))...
                * ((a - A(probe_idx)) .^ 2);

x = 1e160:1e159:1e162;
y = x;
y1 = x;
y2 = x;
for j = 1:length(x)
    y(j) = f_a_plot(x(j));
    y1(j) = f_a_approx(x(j));
end
figure;
plot(x, y, 'r');
hold on;
%plot(x, y1);

%% B quad approx
probe_idx = 18;

Q = langmuir_func(A, B, C);
f_b_plot = @(b) nmf_alpha_beta_divergence(I(probe_idx, :), (A(probe_idx) * C) ./ (1 + b * C), alpha, beta) + 0.5*alpha_B/2 * b^2;

%rpk = (A(array_idx) * C ./ (I(array_idx, :)));
%f_b_approx = @(b) sum(rpk ./ (1 + B_old(array_idx) * C) - (rpk .* C ./ (1 + B_old(array_idx) * C) .^ 2) * (b - B_old(array_idx)) + ...
%    ((rpk .* C .^ 2) ./ (1 + B_old(array_idx) * C) .^ 3) * (b - B_old(array_idx)) .^ 2 + ...
%    log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)) ...
%    - 0.5 * (C .^ 2 ./ (1 + C * B_old(array_idx)) .^ 2) * (b - B_old(array_idx)) .^ 2);
F1 = (A(probe_idx) * C) ./ (1 + B(probe_idx) * C);
f_b_approx = @(b) nmf_alpha_beta_divergence(I(probe_idx, :), (A(probe_idx) * C) ./ (1 + b * C), alpha, beta) + 0.5*alpha_B/2 * b^2 + ...
    ((1 ./ (A(probe_idx) + eps)) .* sum(((F1 + eps) .^ (beta + 1)) .* ((I(probe_idx, :) + eps) .^ alpha - (F1 + eps) .^ alpha), 2) + alpha_B * b) * (b - B(probe_idx)) + ...
                0.5 * (sum(bsxfun(@rdivide, C .^ 2, (1 + B(probe_idx) * C) .^ 2) .* ((F1 + eps) .^ beta) .* ((beta + 1) * ((I(probe_idx, :) + eps) .^ alpha) - ...
                (alpha + beta + 1) * ((F1 + eps) .^ alpha)), 2) + alpha_B) ...
                * ((b - B(probe_idx)) .^ 2);
            
f_b_fd = @(b) ((1 ./ (A(probe_idx) + eps)) .* sum((((A(probe_idx) * C) ./ (1 + b * C) + eps) .^ (beta + 1)) .* ((I(probe_idx, :) + eps) .^ alpha - ((A(probe_idx) * C) ./ (1 + b * C) + eps) .^ alpha), 2) + alpha_B * b);
            
f_b_sd = @(b) 0.5 * (sum(bsxfun(@rdivide, C .^ 2, (1 + b * C) .^ 2) .* (((A(probe_idx) * C) ./ (1 + b * C) + eps) .^ beta) .* ((beta + 1) * ((I(probe_idx, :) + eps) .^ alpha) - ...
                (alpha + beta + 1) * (((A(probe_idx) * C) ./ (1 + b * C) + eps) .^ alpha)), 2) + alpha_B);

x = 0.347:0.00001:0.353;
%x = 0.2:0.001:0.5;
y = x;
y1 = x;
y2 = x;
y3 = x;
for j = 1:length(x)
    y(j) = f_b_plot(x(j));
    y1(j) = f_b_approx(x(j));
    y2(j) = f_b_sd(x(j));
    y3(j) = f_b_fd(x(j));
end
figure;
plot(x, y, 'r');
hold on;
plot(x, y1);
%plot(x, y3, 'k');
%hold on;
%grid on;
%plot([B(probe_idx) B(probe_idx)], [min(y3) max(y3)], '--');