% mad_A = zeros(length(alpha_range), length(beta_range));
% outliers_A = zeros(length(alpha_range), length(beta_range));
% mad_B = zeros(length(alpha_range), length(beta_range));
% outliers_B = zeros(length(alpha_range), length(beta_range));
% C_loo = zeros(length(alpha_range), length(beta_range));
% test_error = zeros(length(alpha_range), length(beta_range));
% overfitting = zeros(length(alpha_range), length(beta_range));

maxIterCnt = 1000;
eps = 1e-5;

for i = 1:length(alpha_range)
    for j = 1:length(beta_range)
        alpha = alpha_range(i);
        beta = beta_range(j);
        fprintf('Alpha = %f, Beta = %f\n', alpha, beta);
        
        if (alpha >= -beta || beta < 0)
            continue;
        end
        
        Abestsliced = A_sliced_all{i, j};
        Bbestsliced = B_sliced_all{i, j};
        Cbest = C_all{i, j};
        
        if (isempty(Abestsliced))
            continue;
        end
        
        [~, test_error(i, j)] = nonlinear_find_concentrations_witherror(inten_test_sliced, Abestsliced, Bbestsliced, ...
            @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
            alpha, beta, maxIterCnt, eps, 10^reg_best(i, j), 1));
        
        overfitting(i, j) = test_error(i, j) - train_error(i, j);
        
        C_loo(i, j) = nonlinear_get_LOO_error(inten_train_sliced, Abestsliced, Bbestsliced, Cbest, reg_best(i, j), ...
            alpha, beta, maxIterCnt, eps);
        
        
        fprintf('Factorizing partitions of I...');
        [mad_A(i, j), outliers_A(i, j), mad_B(i, j), outliers_B(i, j)] = ...
            get_quality_cv_frankenstein(inten_train, inten_train_sliced, inten_full_idx, inten_test_sliced, ...
            partition_arrays, [], [], [], [], ...
            alpha, beta, maxIterCnt, eps, reg_best(i, j));
        fprintf(' Done\n');
        
        %fprintf('Factorizing I...');
        %[A, C, Avect, A_sliced, converged_cnt_full(i, j), norm_problems_cnt(i, j)] = calibrate_model_parallel(inten, inten_sliced, ...
        %    inten_full_idx, @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6));
        %fprintf(' Done\n');
        
        %fprintf('Factorizing I_test with fixed A...');
        %[C_test converged_cnt_fixedA(i, j)] = find_concentrations(inten_test_sliced, A_sliced, ...
        %    @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, maxIterCnt, 1e-6));
        %fprintf(' Done\n');
        
        save('cv_frankenshtein_res_new_crit.mat', 'test_error', 'overfitting', 'mad_A', ...
            'outliers_A', 'mad_B', 'outliers_B', 'C_loo');
    end
end