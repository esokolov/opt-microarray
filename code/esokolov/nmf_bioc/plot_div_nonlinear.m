%% plot divergence(B)
%f1 = @(b) sum(log(I(16, :)) -1 - log(A(16) * C ./ (1 + b * C)) + A(16) * C ./ (I(16, :) .* (1 + b * C)));
array_idx = 1;
f_b_plot = @(b) nmf_alpha_beta_divergence(I(array_idx, :), A(array_idx) * C ./ (1 + b * C), alpha, beta);
x = 0:0.000001:0.001;
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
f_b_approx = @(b) sum(rpk ./ (1 + B_old(array_idx) * C) - (rpk .* C ./ (1 + B_old(array_idx) * C) .^ 2) * (b - B_old(array_idx)) + ...
    ((rpk .* C .^ 2) ./ (1 + B_old(array_idx) * C) .^ 3) * (b - B_old(array_idx)) .^ 2 + ...
    log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)) ...
    - 0.5 * (C .^ 2 ./ (1 + C * B_old(array_idx)) .^ 2) * (b - B_old(array_idx)) .^ 2);

f_b_lin = @(b) sum(rpk ./ (1 + b * C) + log(1 + C * B_old(array_idx)) + (C ./ (1 + C * B_old(array_idx))) * (b - B_old(array_idx)));

x = 0:0.00001:0.001;
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
plot(x, y2, 'c');