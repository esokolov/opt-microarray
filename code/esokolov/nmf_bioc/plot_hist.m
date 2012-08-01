alpha_A_range = [0 10 .^ (-2:2)];
alpha_B_range = [0 10 .^ (-2:2)];
alpha_C_range = 10 .^ (-10:2:2);

qual_hist_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);
C_max_hist_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);
B_max_hist_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);
A_max_hist_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);
corr_B_hist_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);
test_err_all = zeros([length(alpha_A_range), length(alpha_B_range), length(alpha_C_range), 5000]);

for alpha_A_idx = 1:length(alpha_A_range)
    alpha_A = alpha_A_range(alpha_A_idx);
    
    qual_hist_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    C_max_hist_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    B_max_hist_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    A_max_hist_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    corr_B_hist_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    test_err_curr = zeros(length(alpha_B_range) * length(alpha_C_range), 5000);
    parfor i = 0:(length(alpha_B_range) * length(alpha_C_range) - 1)
        alpha_B_idx = mod(i, length(alpha_B_range)) + 1;
        alpha_C_idx = (i - mod(i, length(alpha_B_range))) / length(alpha_B_range) + 1;
        
        alpha_B = alpha_B_range(alpha_B_idx);
        alpha_C = alpha_C_range(alpha_C_idx);
        
        fprintf('%e %e %e\n', alpha_A, alpha_B, alpha_C);

        [A B C isConverged qual_hist_curr(i + 1, :) C_max_hist_curr(i + 1, :) ...
            B_max_hist_curr(i + 1, :) A_max_hist_curr(i + 1, :) corr_B_hist_curr(i + 1, :) test_err_curr(i + 1, :)] = ...
            nonlinear_alpha_beta_linesearch(I(:, 1:500), -0.5, -0.5, 5000, 1e-12, alpha_A, alpha_B, alpha_C, 0, I(:, 501:end));
            %nonlinear_alpha_beta_linesearch(I, -0.5, -0.5, 5000, 1e-12, alpha_A, alpha_B, alpha_C, 0);
    end
    
    for i = 0:(length(alpha_B_range) * length(alpha_C_range) - 1)
        alpha_B_idx = mod(i, length(alpha_B_range)) + 1;
        alpha_C_idx = (i - mod(i, length(alpha_B_range))) / length(alpha_B_range) + 1;
        
        qual_hist_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = qual_hist_curr(i + 1, :);
        C_max_hist_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = C_max_hist_curr(i + 1, :);
        B_max_hist_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = B_max_hist_curr(i + 1, :);
        A_max_hist_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = A_max_hist_curr(i + 1, :);
        corr_B_hist_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = corr_B_hist_curr(i + 1, :);
        test_err_all(alpha_A_idx, alpha_B_idx, alpha_C_idx, :) = test_err_curr(i + 1, :);
    end
    
    save('plot_hist_data.mat', 'qual_hist_all', 'C_max_hist_all', 'B_max_hist_all', 'A_max_hist_all', 'corr_B_hist_all', 'test_err_all');
end

qual_hist_all = shiftdim(qual_hist_all, 3);
C_max_hist_all = shiftdim(C_max_hist_all, 3);
B_max_hist_all = shiftdim(B_max_hist_all, 3);
A_max_hist_all = shiftdim(A_max_hist_all, 3);
corr_B_hist_all = shiftdim(corr_B_hist_all, 3);
test_err_all = shiftdim(test_err_all, 3);

% for alpha_A_idx = 1:length(alpha_A_range)
%    for alpha_B_idx = 1:length(alpha_B_range)
%        for alpha_C_idx = 1:length(alpha_C_range)
%             alpha_A = alpha_A_range(alpha_A_idx);
%             alpha_B = alpha_B_range(alpha_B_idx);
%             alpha_C = alpha_C_range(alpha_C_idx);
%             
%             figure;
%             subplot(4, 1, 1);
%             plot((A_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'b', 'LineWidth', 2);
%             legend('A_{mean}', 'Location', 'NorthWest');
%             title(['alpha_A = ' num2str(alpha_A) '; alpha_B = ' num2str(alpha_B) '; alpha_C = ' num2str(alpha_C)]);
%             grid on;
%             subplot(4, 1, 2);
%             plot((B_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'r', 'LineWidth', 2);
%             legend('B_{mean}', 'Location', 'NorthWest');
%             grid on;
%             subplot(4, 1, 3);
%             plot((C_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'k', 'LineWidth', 2);
%             legend('C_{mean}', 'Location', 'NorthWest');
%             grid on;
%             subplot(4, 1, 4);
%             plot((qual_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
%             legend('Q', 'Location', 'NorthWest');
%             grid on;
% 
%             set(gcf, 'Position', [0 0 1000 700]);
%             %saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.fig']);
%             saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.png'], 'png');
%             %saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.eps']);
%        end
%    end
% end

for alpha_A_idx = 1:length(alpha_A_range)
       parfor i = 0:(length(alpha_B_range) * length(alpha_C_range) - 1)
            alpha_B_idx = mod(i, length(alpha_B_range)) + 1;
            alpha_C_idx = (i - mod(i, length(alpha_B_range))) / length(alpha_B_range) + 1;
            alpha_A = alpha_A_range(alpha_A_idx);
            alpha_B = alpha_B_range(alpha_B_idx);
            alpha_C = alpha_C_range(alpha_C_idx);
            
            figure;
            subplot(6, 1, 1);
            plot((A_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'b', 'LineWidth', 2);
            legend('A_{mean}', 'Location', 'NorthWest');
            title(['alpha_A = ' num2str(alpha_A) '; alpha_B = ' num2str(alpha_B) '; alpha_C = ' num2str(alpha_C)]);
            grid on;
            subplot(6, 1, 2);
            plot((B_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'r', 'LineWidth', 2);
            legend('B_{mean}', 'Location', 'NorthWest');
            grid on;
            subplot(6, 1, 3);
            plot((C_max_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'k', 'LineWidth', 2);
            legend('C_{mean}', 'Location', 'NorthWest');
            grid on;
            subplot(6, 1, 4);
            plot((qual_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
            legend('Q', 'Location', 'NorthWest');
            grid on;
            subplot(6, 1, 5);
            plot((test_err_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'Color', [148/256 0/256 211/256], 'LineWidth', 2);
            legend('Q_{test}', 'Location', 'NorthWest');
            grid on;
            subplot(6, 1, 6);
            plot((corr_B_hist_all(:, alpha_A_idx, alpha_B_idx, alpha_C_idx)), 'Color', [46/256 139/256 87/256], 'LineWidth', 2);
            legend('Corr(B, I)', 'Location', 'NorthWest');
            grid on;

            set(gcf, 'Position', [0 0 1000 900]);
            saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.fig'], 'fig');
            saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.png'], 'png');
            saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.eps'], 'eps');
       end
end