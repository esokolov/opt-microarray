load('probeset_6519_all.mat');

alpha_reg_range = 10 .^ (-10:0);


qual_hist_all = zeros([length(alpha_reg_range), 5000]);
C_max_hist_all = zeros([length(alpha_reg_range),5000]);
B_max_hist_all = zeros([length(alpha_reg_range), 5000]);
A_max_hist_all = zeros([length(alpha_reg_range), 5000]);
corr_B_hist_all = zeros([length(alpha_reg_range), 5000]);
test_err_all = zeros([length(alpha_reg_range), 5000]);

parfor alpha_reg_idx = 1:length(alpha_reg_range)
    alpha_reg = alpha_reg_range(alpha_reg_idx);
                
    fprintf('%e\n', alpha_reg);

    [A B C isConverged qual_hist_all(alpha_reg_idx, :) C_max_hist_all(alpha_reg_idx, :) ...
        B_max_hist_all(alpha_reg_idx, :) A_max_hist_all(alpha_reg_idx, :) ...
        corr_B_hist_all(alpha_reg_idx, :) test_err_all(alpha_reg_idx, :)] = ...
        nonlinear_alpha_beta_reg_derivative_4pics(I(:, 1:500), -0.5, -0.5, 5000, 1e-6, alpha_reg, 0, I(:, 501:end));
        %nonlinear_alpha_beta_linesearch(I, -0.5, -0.5, 5000, 1e-12, alpha_A, alpha_B, alpha_C, 0);
    fprintf('%d completed\n', alpha_reg_idx);
end

save('plot_hist_reg_data_6519.mat', 'qual_hist_all', 'C_max_hist_all', 'B_max_hist_all', 'A_max_hist_all', 'corr_B_hist_all', 'test_err_all');

parfor alpha_reg_idx = 1:length(alpha_reg_range)
        alpha_reg = alpha_reg_range(alpha_reg_idx);

        figure;
        subplot(6, 1, 1);
        plot((A_max_hist_all(alpha_reg_idx, :)), 'b', 'LineWidth', 2);
        legend('A_{mean}', 'Location', 'NorthWest');
        title(['alpha_{reg} = ' num2str(alpha_reg)]);
        grid on;
        subplot(6, 1, 2);
        plot((B_max_hist_all(alpha_reg_idx, :)), 'r', 'LineWidth', 2);
        legend('B_{mean}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 3);
        plot((C_max_hist_all(alpha_reg_idx, :)), 'k', 'LineWidth', 2);
        legend('C_{mean}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 4);
        plot((qual_hist_all(alpha_reg_idx, :)), 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
        legend('Q', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 5);
        plot((test_err_all(alpha_reg_idx, :)), 'Color', [148/256 0/256 211/256], 'LineWidth', 2);
        legend('Q_{test}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 6);
        plot((corr_B_hist_all(alpha_reg_idx, :)), 'Color', [46/256 139/256 87/256], 'LineWidth', 2);
        legend('Corr(B, I)', 'Location', 'NorthWest');
        grid on;

        set(gcf, 'Position', [0 0 1000 900]);
        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.fig'], 'fig');
        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.png'], 'png');
        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.eps'], 'eps');
end