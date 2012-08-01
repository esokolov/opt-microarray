alpha_range = -4:0.25:4;
beta_range = -4:0.25:4;

goodness_of_fit_huber = zeros(length(alpha_range), length(beta_range), length(I_10));
corr_B = zeros(length(alpha_range), length(beta_range), length(I_10));

for I_idx = 1:length(I_10)
    for i = 1:length(alpha_range)
        parfor j = 1:length(beta_range)
            alpha = alpha_range(i);
            beta = beta_range(j);
            fprintf('%f %f\n', alpha, beta);
            %[A B C] = nonlinear_alpha_beta_reg_derivative(I_10{I_idx}, alpha, beta, 500, 1e-12, 1e-4, 0);
			[A B C] = nonlinear_alpha_beta_linesearch(I_10{I_idx}, alpha, beta, 500, 1e-12, 0, 0.01, 1e-10, 0);
            corr_B(i, j, I_idx) = corr(B, quantile(I_10{I_idx}', 0.9)', 'type', 'Spearman');
            goodness_of_fit_huber(i, j, I_idx) = sum(sum(huber_func(I_10{I_idx} - langmuir_func(A, B, C))));
        end
    end
    save('corr_and_huber_10_quad', 'goodness_of_fit_huber', 'corr_B');
end

save('corr_and_huber_10_quad', 'goodness_of_fit_huber', 'corr_B');
