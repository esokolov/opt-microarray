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

%fprintf('Generating partitions...');
%[partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
%   generate_partitions(inten_train, inten_train_sliced, inten_full_idx);
%fprintf(' Done\n');
load('partition_arrays.mat');

maxIterCnt = 1000;
eps = 1e-5;

[A_train, C_train, ~, A_train_sliced] = calibrate_model_parallel(inten_train, inten_train_sliced, inten_full_idx(1:50), @(I) nmf_rma(I));
[A_test_1, C_test_1, ~, A_test_sliced_1] = calibrate_model_parallel(inten_test(1:1568, :), inten_test_sliced, inten_full_idx(1:50), @(I) nmf_rma(I));
[C_test_2] = find_concentrations(inten_test_sliced, A_train_sliced, @(I, A) nmf_rma_fixedA(I, A));

[A_train_LO, C_train_LO, ~, A_train_sliced_LO] = calibrate_model_parallel(inten_train, inten_train_sliced, ...
    inten_full_idx(1:50), @(I) nmf_rma_LO_reg(I));
[C_test_LO] = find_concentrations(inten_test_sliced, ...
    A_train_sliced_LO, @(I, A) nmf_rma_LO_reg_fixedA(I, A));

LOO_error_rma = nmf_rma_get_LOO_error(inten_train_sliced, A_train_sliced_LO, C_train_LO);

train_error_rma = quality_functional_linear(inten_train_sliced, A_train_sliced, C_train);
test_error_rma_1 = quality_functional_linear(inten_test_sliced, A_test_sliced_1, C_test_1);
test_error_rma_2 = quality_functional_linear(inten_test_sliced, A_train_sliced, C_test_2);
test_error_rma_LO = quality_functional_linear(inten_test_sliced, A_train_sliced_LO, C_test_LO);

fprintf('Factorizing partitions of I...');
[mad_A_rma, outliers_A_rma] = ...
    get_quality_cv_rma(inten_train, inten_train_sliced, inten_full_idx, partition_arrays);
fprintf(' Done\n');


save('cv_rma_new_res.mat', 'mad_A_rma', 'outliers_A_rma', ...
    'train_error_rma', 'test_error_rma_1', 'test_error_rma_2', 'test_error_rma_LO', 'LOO_error_rma');
