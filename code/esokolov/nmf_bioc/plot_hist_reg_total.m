%load('probeset_6519_all.mat');

alpha_reg_range = [0 10 .^ (-7:2)];
maxIterCnt = 1000;
alpha = 3;
beta = -1;
probesets = 1:100;

qual_hist_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);
C_max_hist_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);
B_max_hist_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);
A_max_hist_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);
corr_B_hist_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);
test_err_all = zeros([length(alpha_reg_range)*length(probesets), maxIterCnt]);

parfor ind = 1:(length(alpha_reg_range)*length(probesets))
    probeset_idx = probesets(mod(ind-1,length(probesets))+1);
    alpha_reg_idx = round((ind-probeset_idx)/length(probesets)) + 1;
    alpha_reg = alpha_reg_range(alpha_reg_idx);
    %for alpha_reg_idx = 1:length(alpha_reg_range)
    %    alpha_reg = alpha_reg_range(alpha_reg_idx);
    
    %    fprintf('%e\n', alpha_reg);
    %    parfor probeset_idx = 1:length(probesets)
    tic;
    I = inten_full_sliced{probeset_idx};
    [A B C isConverged qual_hist_all(ind, :) C_max_hist_all(ind, :) ...
        B_max_hist_all(ind, :) A_max_hist_all(ind, :) ...
        corr_B_hist_all(ind, :) test_err_all(ind, :)] = ...
        nonlinear_alpha_beta_reg_derivative_4pics(I(:, 1:200), alpha, beta, maxIterCnt, 1e-6, alpha_reg, 0, I(:, 801:end));
    fprintf('%e - %d; ', alpha_reg,probeset_idx);
    toc;
    %    end
    %nonlinear_alpha_beta_linesearch(I, -0.5, -0.5, 5000, 1e-12, alpha_A, alpha_B, alpha_C, 0);
    %    fprintf('%d completed\n', alpha_reg_idx);
    %end
end



save('plot_hist_reg_data_200arrays_100probesets.mat', 'qual_hist_all', 'C_max_hist_all', 'B_max_hist_all', 'A_max_hist_all', 'corr_B_hist_all', 'test_err_all');

qual_hist_total = zeros([length(alpha_reg_range), maxIterCnt]);
C_max_hist_total = zeros([length(alpha_reg_range), maxIterCnt]);
B_max_hist_total = zeros([length(alpha_reg_range), maxIterCnt]);
A_max_hist_total = zeros([length(alpha_reg_range), maxIterCnt]);
corr_B_hist_total = zeros([length(alpha_reg_range), maxIterCnt]);
test_err_total = zeros([length(alpha_reg_range), maxIterCnt]);

for ind = 1:length(alpha_reg_range)
    tosum_ind = length(probesets)*repmat(ind-1,1,length(probesets)) + [1:length(probesets)];
    qual_hist_total(ind,:) = sum(qual_hist_all(tosum_ind,:),1);
    C_max_hist_total(ind,:) = mean(C_max_hist_all(tosum_ind,:),1);
    B_max_hist_total(ind,:) = mean(B_max_hist_all(tosum_ind,:),1);
    A_max_hist_total(ind,:) = mean(A_max_hist_all(tosum_ind,:),1);
    corr_B_hist_total(ind,:) = mean(corr_B_hist_all(tosum_ind,:),1);
    test_err_total(ind,:) = sum(test_err_all(tosum_ind,:),1);
end


parfor alpha_reg_idx = 1:length(alpha_reg_range)
        alpha_reg = alpha_reg_range(alpha_reg_idx);

        figure;
        subplot(6, 1, 1);
        plot((A_max_hist_total(alpha_reg_idx, :)), 'b', 'LineWidth', 2);
        legend('A_{mean}', 'Location', 'NorthWest');
        title(['alpha_{reg} = ' num2str(alpha_reg)]);
        grid on;
        subplot(6, 1, 2);
        plot((B_max_hist_total(alpha_reg_idx, :)), 'r', 'LineWidth', 2);
        legend('B_{mean}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 3);
        plot((C_max_hist_total(alpha_reg_idx, :)), 'k', 'LineWidth', 2);
        legend('C_{mean}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 4);
        plot((qual_hist_total(alpha_reg_idx, :)), 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
        legend('Q', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 5);
        plot((test_err_total(alpha_reg_idx, :)), 'Color', [148/256 0/256 211/256], 'LineWidth', 2);
        legend('Q_{test}', 'Location', 'NorthWest');
        grid on;
        subplot(6, 1, 6);
        plot((corr_B_hist_total(alpha_reg_idx, :)), 'Color', [46/256 139/256 87/256], 'LineWidth', 2);
        legend('Corr(B, I)', 'Location', 'NorthWest');
        grid on;

        set(gcf, 'Position', [0 0 1000 900]);
        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.fig'], 'fig');
        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.png'], 'png');
%        saveas(gcf, ['reg_pics/reg_' int2str(alpha_reg_idx) '.eps'], 'eps');
end