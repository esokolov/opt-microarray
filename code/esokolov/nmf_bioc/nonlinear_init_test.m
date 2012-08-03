[A_nmf, B_nmf, C_nmf, ~, qual_hist_nmf] = nonlinear_alpha_beta_reg_derivative(I, -0.5, -0.5, 500, 1e-6, 1e-5, 0);
qual_nmf = nmf_alpha_beta_divergence(I, langmuir_func(A_nmf, B_nmf, C_nmf), -0.5, -0.5);

%%
quals = zeros(10, 1);
qual_hists = cell(10, 1);
for i = 1:10
    [A, B, C, ~, qual_hists{i}] = nonlinear_alpha_beta_reg_derivative_rand_init(I, -0.5, -0.5, 500, 1e-6, 1e-5, 0);
    quals(i) = nmf_alpha_beta_divergence(I, langmuir_func(A, B, C), -0.5, -0.5);
    fprintf('%d %e\n', i, quals(i));
end

%%
cc = hsv(11);

figure;
plot(1:length(qual_hist_nmf)-1, qual_hist_nmf(1:end-1), 'LineWidth', 3, 'color', cc(1, :));
hold on;
for i = 1:10
   plot(1:length(qual_hists{i})-1, qual_hists{i}(1:end-1), 'LineWidth', 2, 'color', cc(1+i, :));
end
