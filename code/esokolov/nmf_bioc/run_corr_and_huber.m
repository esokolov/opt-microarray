alpha_range = -10:0.25:10;
beta_range = -10:0.25:10;

%goodness_of_fit_huber = zeros(length(alpha_range), length(beta_range));
%corr_B = zeros(length(alpha_range), length(beta_range));

for i = 41:length(alpha_range)
    parfor j = 1:length(beta_range)
        alpha = alpha_range(i);
        beta = beta_range(j);
        fprintf('%f %f\n', alpha, beta);
        [A B C] = nonlinear_alpha_beta_reg_derivative(I, alpha, beta, 500, 1e-12, 1e-4, 0);
        corr_B(i, j) = corr(B, quantile(I', 0.9)', 'type', 'Spearman');
        goodness_of_fit_huber(i, j) = sum(sum(huber_func(I - langmuir_func(A, B, C))));
    end
end

save('corr_and_huber', 'goodness_of_fit_huber', 'corr_B');