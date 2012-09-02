js = [15,15,15,15,14,17,14];
is = [4, 5, 6, 8, 7, 4, 4];

for k=1:length(js)
    i = is(k);
    j = js(k);
    alpha = alpha_range(i);
    beta = beta_range(j);
    fprintf('Alpha = %f, Beta = %f\n', alpha, beta);
    
    if (alpha < -beta)
        continue;
    end
    
    tic;
    [Abest, Bbest, Cbest, Cbestcontrol, Abestsliced, Bbestsliced, reg_best(i, j), test_error(i, j)] = ...
        nonlinear_calibrate_frankenstein_model_GS(inten_train, inten_train_sliced, inten_test_sliced, inten_full_idx, ...
        alpha, beta, maxIterCnt, eps);
    
    for gene_idx = 1:length(inten_train_sliced)
        I = inten_train_sliced{gene_idx};
        error = (langmuir_func(Abestsliced{gene_idx}, Bbestsliced{gene_idx}, Cbest(gene_idx, :)) - I) ./ ...
            I .* repmat(Cbest(gene_idx, :), size(I, 1), 1);
        bound = quantile(error(:), 0.75);
        W = ones(size(I));
        W = W .* (error <= bound);
        train_error(i, j) = train_error(i, j) + ...
            sum(sum(I .* abs(langmuir_func(Abestsliced{gene_idx}, Bbestsliced{gene_idx}, Cbest(gene_idx, :)) - I) .* W)) / train_size;
    end
    
    overfitting(i, j) = test_error(i, j) - train_error(i, j);
    
    C_loo(i, j) = nonlinear_get_LOO_error(inten_train_sliced, Abestsliced, Bbestsliced, Cbest, reg_best(i, j), ...
        alpha, beta, maxIterCnt, eps);
    
    
    fprintf('Factorizing partitions of I...');
    [mad_A(i, j), outliers_A(i, j), mad_B(i, j), outliers_B(i, j)] = ...
        get_quality_cv_frankenstein(inten_train, inten_train_sliced, inten_full_idx, inten_test_sliced, ...
        partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl, ...
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
    
    A_all{i, j} = Abest;
    B_all{i, j} = Bbest;
    C_all{i, j} = Cbest;
    
    A_sliced_all{i, j} = Abestsliced;
    B_sliced_all{i, j} = Bbestsliced;
    
    save('cv_frankenstein_res_bigger.mat', 'mad_A', 'outliers_A', 'mad_B', 'outliers_B', 'C_loo', ...
        'train_error', 'test_error', 'overfitting', 'A_all', 'B_all', 'C_all', 'A_sliced_all', 'B_sliced_all', 'reg_best');
    toc;
end