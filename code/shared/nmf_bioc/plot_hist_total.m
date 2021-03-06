%alpha_A_range = [0 10 .^ (-2:2)];
alpha_B_range = [0 10 .^ (-3:3)];
alpha_C_range = [0 10 .^ (-6:2:2)];
maxIterCnt = 1000;
alpha = 3;
beta = -1;
probesets = 1:50;

qual_hist_all   = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);
C_max_hist_all  = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);
B_max_hist_all  = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);
A_max_hist_all  = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);
corr_B_hist_all = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);
test_err_all    = zeros([length(alpha_B_range), length(alpha_C_range), maxIterCnt]);

qual_hist_curr   = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);
C_max_hist_curr  = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);
B_max_hist_curr  = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);
A_max_hist_curr  = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);
corr_B_hist_curr = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);
test_err_curr    = zeros(length(probesets) * length(alpha_C_range), maxIterCnt);

for alpha_B_idx = 1:length(alpha_B_range)
    alpha_B = alpha_B_range(alpha_B_idx);
    
    fprintf('alpha_B:%e \n', alpha_B);
        
    parfor ind = 1:(length(alpha_C_range)*length(probesets))
        tic;
        probeset_idx = probesets(mod(ind-1,length(probesets))+1);
        alpha_C_idx = round((ind-(mod(ind-1,length(probesets))+1))/length(probesets)) + 1;
        alpha_C = alpha_C_range(alpha_C_idx);
        I = inten_full_sliced{probeset_idx};                
        
        [A B C isConverged qual_hist_curr(ind, :) C_max_hist_curr(ind, :) ...
            B_max_hist_curr(ind, :) A_max_hist_curr(ind, :) corr_B_hist_curr(ind, :) test_err_curr(ind, :)] = ...
            nonlinear_alpha_beta_linesearch_4pics(I(:, 1:100), alpha, beta, maxIterCnt, 1e-6, 0, alpha_B, alpha_C, 0, I(:, 901:end));
        
        fprintf('probeset %d, alpha_C=%e. ', probeset_idx, alpha_C); 
        toc;
    end
    
    for ind = 1:length(alpha_C_range)
        tosum_ind = length(probesets)*repmat(ind-1,1,length(probesets)) + [1:length(probesets)];
        qual_hist_all(alpha_B_idx, ind, :)   = sum(qual_hist_curr(tosum_ind,:), 1);
        C_max_hist_all(alpha_B_idx, ind, :)  = mean(C_max_hist_curr(tosum_ind,:), 1);
        B_max_hist_all(alpha_B_idx, ind, :)  = mean(B_max_hist_curr(tosum_ind,:),1);
        A_max_hist_all(alpha_B_idx, ind, :)  = mean(A_max_hist_curr(tosum_ind,:), 1);
        corr_B_hist_all(alpha_B_idx, ind, :) = mean(corr_B_hist_curr(tosum_ind,:), 1);
        test_err_all(alpha_B_idx, ind, :)    = sum(test_err_curr(tosum_ind,:), 1);
    end
    
    save('plot_hist_data_200arrays_100probesets.mat', 'qual_hist_all', 'C_max_hist_all', 'B_max_hist_all', 'A_max_hist_all', 'corr_B_hist_all', 'test_err_all');
end

parfor i = 0:(length(alpha_B_range) * length(alpha_C_range) - 1)
    alpha_B_idx = mod(i, length(alpha_B_range)) + 1;
    alpha_C_idx = (i - mod(i, length(alpha_B_range))) / length(alpha_B_range) + 1;
    alpha_B = alpha_B_range(alpha_B_idx);
    alpha_C = alpha_C_range(alpha_C_idx);
    
    figure;
    subplot(6, 1, 1);
    tmp = (A_max_hist_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'b', 'LineWidth', 2);
    legend('A_{mean}', 'Location', 'NorthWest');
    title(['alpha_B = ' num2str(alpha_B) '; alpha_C = ' num2str(alpha_C)]);
    grid on;
    
    subplot(6, 1, 2);
    tmp = (B_max_hist_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'r', 'LineWidth', 2);
    legend('B_{mean}', 'Location', 'NorthWest');
    grid on;
    
    subplot(6, 1, 3);
    tmp = (C_max_hist_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'k', 'LineWidth', 2);
    legend('C_{mean}', 'Location', 'NorthWest');
    grid on;
    
    subplot(6, 1, 4);
    tmp = (qual_hist_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
    legend('Q', 'Location', 'NorthWest');
    grid on;
    
    subplot(6, 1, 5);
    tmp = (test_err_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'Color', [148/256 0/256 211/256], 'LineWidth', 2);
    legend('Q_{test}', 'Location', 'NorthWest');
    grid on;
    
    subplot(6, 1, 6);
    tmp = (corr_B_hist_all(alpha_B_idx, alpha_C_idx, :)); 
    plot(tmp(:), 'Color', [46/256 139/256 87/256], 'LineWidth', 2);
    legend('Corr(B, I)', 'Location', 'NorthWest');
    grid on;
    
    set(gcf, 'Position', [0 0 1000 900]);
    saveas(gcf, ['reg_pics/reg_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.fig'], 'fig');
    saveas(gcf, ['reg_pics/reg_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.png'], 'png');
    %saveas(gcf, ['reg_pics/reg_' int2str(alpha_A_idx) '_' int2str(alpha_B_idx) '_' int2str(alpha_C_idx) '.eps'], 'eps');
end