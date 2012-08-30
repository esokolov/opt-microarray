alpha_range = -3:0.5:4;
beta_range = -2:0.5:4;
%alpha_range = -2:0.25:2;
%beta_range = -2:0.25:4;

mad_A = zeros(length(alpha_range), length(beta_range));
ouliers_A = zeros(length(alpha_range), length(beta_range));
mad_B = zeros(length(alpha_range), length(beta_range));
ouliers_B = zeros(length(alpha_range), length(beta_range));
C_loo = zeros(length(alpha_range), length(beta_range));
train_error = zeros(length(alpha_range), length(beta_range));
test_error = zeros(length(alpha_range), length(beta_range));
overfitting = zeros(length(alpha_range), length(beta_range));
reg_best = zeros(length(alpha_range), length(beta_range));

A_all = cell(length(alpha_range), length(beta_range));
B_all = cell(length(alpha_range), length(beta_range));
C_all = cell(length(alpha_range), length(beta_range));

%inten = inten_full;
%inten_sliced = inten_full_sliced;

% inten = inten_full(:, [1:100 1001:1002]) + 1;
% inten_sliced = inten_full_sliced;
% for i = 1:length(inten_sliced)
%     inten_sliced{i} = inten_sliced{i}(:, 1:100) + 1;
% end
% 
% inten_test = inten_full(:, 901:end) + 1;
% inten_test_sliced = inten_full_sliced;
% for i = 1:length(inten_test_sliced)
%     inten_test_sliced{i} = inten_test_sliced{i}(:, 901:end) + 1;
% end

train_size = 200;
test_size = 200;

%formsamples;

fprintf('Generating partitions...');
[partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
   generate_partitions(inten_train, inten_train_sliced, inten_full_idx);
fprintf(' Done\n');
%load('alpha_beta_partition.mat');

maxIterCnt = 1000;
eps = 1e-5;

for i = 1:length(alpha_range)
    for j = 1:length(beta_range)
        alpha = alpha_range(i);
        beta = beta_range(j);
        fprintf('Alpha = %f, Beta = %f\n', alpha, beta);
        
        if (alpha < -beta)
            continue;
        end

        fprintf('Factorizing partitions of I...');
        [mad_A, outliers_A, mad_B, outliers_B] = ...
            get_quality_cv_frankenstein(inten_train, inten_train_sliced, inten_full_idx, inten_test_sliced, ...
            partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl, ...
            alpha, beta, maxIterCnt, eps);
        fprintf(' Done\n');
        
        [Abest, Bbest, Cbest, Cbestcontrol, Abestsliced, Bbestsliced, reg_best(i, j), test_error(i, j)] = ...
            nonlinear_calibrate_frankenstein_model(inten_train, inten_train_sliced, inten_test_sliced, inten_full_idx, ...
            alpha, beta, maxIterCnt, eps);
        
        I = inten_train(:, 1:end-2);
        error = (langmuir_func(Abest, Bbest, Cbest) - I) ./ I .* repmat(C, size(I, 1), 1);
        bound = quantile(error(:), 0.75);
        W = ones(size(inten));
        W = W .* (error <= bound);
        train_error(i, j) = sum(sum(I .* abs(langmuir_func(Abest, Bbest, Cbest) - I) .* W)) / train_size;
        
        overfitting(i, j) = test_error(i, j) - train_error(i, j);
        
        C_loo(i, j) = nonlinear_get_LOO_error(inten_train_sliced, Abestsliced, Bbestsliced, Cbest, reg_best(i, j), ...
            alpha, beta, maxIterCnt, eps);
        
        %fprintf('Factorizing I...');
        %[A, C, Avect, A_sliced, converged_cnt_full(i, j), norm_problems_cnt(i, j)] = calibrate_model_parallel(inten, inten_sliced, ...
        %    inten_full_idx, @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6));
        %fprintf(' Done\n');
        
        %fprintf('Factorizing I_test with fixed A...');
        %[C_test converged_cnt_fixedA(i, j)] = find_concentrations(inten_test_sliced, A_sliced, ...
        %    @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, maxIterCnt, 1e-6));
        %fprintf(' Done\n');
        
        save('cv_frankenstein_res.mat', 'mad_A', 'outliers_A', 'mad_B', 'outliers_B', 'C_loo', ...
            'train_error', 'test_error', 'overfitting', 'A_all', 'B_all', 'C_all');
    end
end