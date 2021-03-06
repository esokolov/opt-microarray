%%
matlabpool open Profile1 50

%%
load('/home/affy/model/inten_all.mat');

K = size(inten_full, 2);

inten_test = inten_full(:, 3001:end) + 1;
inten_test_sliced = inten_full_sliced;
for i = 1:length(inten_test_sliced)
    inten_test_sliced{i} = inten_test_sliced{i}(:, 3001:end) + 1;
end

alpha = -0.5;
beta = -0.5;

%%
train_size_range = 100;%:100:3000;

mad_A = zeros(length(train_size_range), 1);
outliers_A = zeros(length(train_size_range), 1);
mad_B = zeros(length(train_size_range), 1);
outliers_B = zeros(length(train_size_range), 1);
mad_C = zeros(length(train_size_range), 1);
outliers_C = zeros(length(train_size_range), 1);
concentration_dists = zeros(length(train_size_range), 1);
goodness_of_fit_native = zeros(length(train_size_range), 1);
goodness_of_fit_huber = zeros(length(train_size_range), 1);
test_error_native = zeros(length(train_size_range), 1);
test_error_huber = zeros(length(train_size_range), 1);
overfitting_native = zeros(length(train_size_range), 1);
overfitting_huber = zeros(length(train_size_range), 1);

fprintf('Generating partitions...');
[~, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
   generate_partitions(inten_full, inten_full_sliced, inten_full_idx);
fprintf(' Done\n');

%%  ��� �������� ��� ������ �������� �������
factorization_func = @(I) nonlinear_alpha_beta_linesearch(I, alpha, beta, 1000, 1e-6, 0, 0, 0, 1);
factorization_func_fixedAB = @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_fixedAB(I_arg, A_arg, B_arg, alpha, beta, 1000, 1e-6, 0, 1);

for train_size_idx = 1:length(train_size_range)
    train_size = train_size_range(train_size_idx);
    fprintf('Training sample size: %d\n', train_size);
    
    inten = inten_full(:, [1:train_size (K - 1):K]) + 1;
    inten_sliced = inten_full_sliced;
    for i = 1:length(inten_sliced)
        inten_sliced{i} = inten_sliced{i}(:, 1:train_size) + 1;
    end
    
   partition_arrays = generate_partitions(inten, inten_sliced, inten_full_idx);
    
    fprintf('Factorizing partitions of I...');
    [mad_A(train_size_idx), outliers_A(train_size_idx), mad_B(train_size_idxi), ...
        outliers_B(train_size_idx), mad_C(train_size_idx), outliers_C(train_size_idx), ...
        arrays_factors, probes_factors] = ...
        get_quality_cv_nonlinear(inten, inten_sliced, inten_full_idx, ...
        factorization_func, partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl);
    fprintf(' Done\n');
    
    conc_gene_dist = zeros(size(probes_factors{3}, 1), 1);
    for k = 1:length(conc_gene_dist)
        conc_gene_dist(k) = Cdist(probes_factors{3}(k, :), probes_factors{10}(k, :));
    end
    concentration_dists(train_size_idx) = mad_my(conc_gene_dist);
    
    fprintf('Factorizing I...');
    [A, B, C, Avect, Bvect, A_sliced, B_sliced] = nonlinear_calibrate_model(inten, inten_sliced, ...
        inten_full_idx, factorization_func, false);
    fprintf(' Done\n');

    fprintf('Factorizing I_test with fixed A...');
    C_test = nonlinear_find_concentrations(inten_test_sliced, A_sliced, B_sliced, factorization_func_fixedAB);
    fprintf(' Done\n');
    
    Qtrain = langmuir_func(A, B, C);
    R = inten(:, 1:end-2) - Qtrain;
    goodness_of_fit_huber(train_size_idx) = sum(sum(huber_func(R)));
    goodness_of_fit_native(train_size_idx) = nmf_alpha_beta_divergence(inten(:, 1:end-2), Qtrain, alpha, beta);
    
    Qtest = langmuir_func(A, B, C_test);
    R = inten_test(:, 1:end-2) - Qtest;
    test_error_huber(train_size_idx) = sum(sum(huber_func(R)));
    test_error_native(train_size_idx) = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), Qtest, alpha, beta);

    
    overfitting_native(train_size_idx) = test_error_native(train_size_idx) / (size(inten_test, 2) - 2) - ...
        goodness_of_fit_native(train_size_idx) / (size(inten, 2) - 2);
    overfitting_huber(train_size_idx) = test_error_huber(train_size_idx) / (size(inten_test, 2) - 2) - ...
        goodness_of_fit_huber(train_size_idx) / (size(inten, 2) - 2);
    
    save('nonlinear_cv_all_arrays.mat', 'mad_A', 'outliers_A', 'mad_B', 'outliers_B', ...
        'mad_C', 'outliers_C', 'concentration_dists', 'goodness_of_fit_native', 'goodness_of_fit_huber', ...
        'test_error_native', 'test_error_huber', 'overfitting_native', 'overfitting_huber');
end

%% �������� ��� RMA
mad_A_rma = zeros(length(train_size_range), 1);
outliers_A_rma = zeros(length(train_size_range), 1);
mad_B_rma = zeros(length(train_size_range), 1);
outliers_B_rma = zeros(length(train_size_range), 1);
mad_C_rma = zeros(length(train_size_range), 1);
outliers_C_rma = zeros(length(train_size_range), 1);
concentration_dists_rma = zeros(length(train_size_range), 1);
goodness_of_fit_huber_rma = zeros(length(train_size_range), 1);

%%
factorization_func = @(I) nmf_rma(I);

for train_size_idx = 1:length(train_size_range)
    train_size = train_size_range(train_size_idx);
    fprintf('Training sample size: %d\n', train_size);
    
    inten = inten_full(:, [1:train_size (K - 1):K]) + 1;
    inten_sliced = inten_full_sliced;
    for i = 1:length(inten_sliced)
        inten_sliced{i} = inten_sliced{i}(:, 1:train_size) + 1;
    end
    
    fprintf('Factorizing partitions of I...');
    [mad_A_rma(train_size_idx), outliers_A_rma(train_size_idx), mad_B_rma(train_size_idxi), ...
        outliers_B_rma(train_size_idx), mad_C_rma(train_size_idx), outliers_C_rma(train_size_idx), ...
        arrays_factors, probes_factors] = ...
        get_quality_cv(inten, inten_sliced, inten_full_idx, ...
        factorization_func, partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl);
    fprintf(' Done\n');
    
    conc_gene_dist = zeros(size(probes_factors{2}, 1), 1);
    for k = 1:length(conc_gene_dist)
        conc_gene_dist(k) = Cdist(probes_factors{2}(k, :), probes_factors{6}(k, :));
    end
    concentration_dists_rma(train_size_idx) = mad_my(conc_gene_dist);
    
    fprintf('Factorizing I...');
    [A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, ...
            inten_full_idx, factorization_func);
    fprintf(' Done\n');
    
    Qtrain = A * C;
    R = inten(:, 1:end-2) - Qtrain;
    goodness_of_fit_huber(train_size_idx) = sum(sum(huber_func(R)));
    
    save('nonlinear_cv_rma_all_arrays.mat', 'mad_A_rma', 'outliers_A_rma', 'mad_B_rma', 'outliers_B_rma', ...
        'mad_C_rma', 'outliers_C_rma', 'concentration_dists_rma', 'goodness_of_fit_huber_rma');
end

    