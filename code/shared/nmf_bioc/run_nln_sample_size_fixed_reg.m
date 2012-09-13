train_sizes = 100:100:2400;
alpha = 2;
beta = 1;
maxIterCnt = 1000;
eps = 1e-5;

alpha_C = -5;

train_error = zeros(length(train_sizes), 1);
test_error = zeros(length(train_sizes), 1);

A_all = cell(length(train_sizes), 1);
A_sliced_all = cell(length(train_sizes), 1);
B_all = cell(length(train_sizes), 1);
B_sliced_all = cell(length(train_sizes), 1);
C_train_all = cell(length(train_sizes), 1);
C_test_all = cell(length(train_sizes), 1);

% train_error_rma = zeros(length(train_sizes), 1);
% test_error_rma_1 = zeros(length(train_sizes), 1);
% test_error_rma_2 = zeros(length(train_sizes), 1);
% test_error_rma_LO = zeros(length(train_sizes), 1);

%arrays_cnt = size(inten_full, 2) - 2;
%validation_arrays = randsample(arrays_cnt, 200);
%test_arrays = randsample(setdiff(1:arrays_cnt, validation_arrays), 500);
%load('test_arrays_idx');

%probes_cnt = inten_full_idx{50}(end);
%inten_validation = inten_full(1:probes_cnt, [validation_arrays' end-1 end]);
%inten_test = inten_full(1:probes_cnt, [test_arrays end-1 end]);
%inten_validation_sliced = inten_full_sliced(1:50);
%inten_test_sliced = inten_full_sliced(1:50);
%for i = 1:50
%    inten_validation_sliced{i} = inten_validation_sliced{i}(:, validation_arrays);
%    inten_test_sliced{i} = inten_test_sliced{i}(:, test_arrays);
%end

%inten_full_idx_curr = inten_full_idx(1:50);

% 
% inten_full_curr = inten_full(1:probes_cnt, [setdiff(1:arrays_cnt, union(validation_arrays, test_arrays)) end-1 end]);
% inten_full_sliced_curr = inten_full_sliced(1:50);
% for i = 1:50
%     inten_full_sliced_curr{i} = inten_full_sliced_curr{i}(:, setdiff(1:arrays_cnt, union(validation_arrays, test_arrays)));
% end
% 
% inten_train_big = zeros(probes_cnt, 1000+2);
% inten_train_big_sliced = inten_full_sliced_curr;
% arrays_train = cell(50, 1);
% for j = 1:50
%    arrays_train{j} = pick_arrays(inten_full_sliced_curr{j}, 1000);
%    inten_train_big_sliced{j} = inten_full_sliced_curr{j}(:, arrays_train{j});
%    inten_train_big(inten_full_idx{j}, :) = inten_full_curr(inten_full_idx{j}, [arrays_train{j} end-1 end]);
%end

for i = 1:length(train_sizes)
    train_size = train_sizes(i);
    
    fprintf('%d\n', train_size);
    
    inten_train_tmp = inten_train(:,[1:train_size end-1 end]);
    inten_train_sliced_tmp = inten_train_sliced;
    for j = 1:50
        inten_train_sliced_tmp{j} = inten_train_sliced{j}(:, 1:train_size);
    end    
    
    [A, B, C_train, ~, ~, A_sliced, B_sliced, train_error(i)] = nonlinear_calibrate_model_witherror(inten_train_tmp, ...
        inten_train_sliced_tmp, inten_full_idx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_C, 1));
    
    [C_test, test_error(i)] = nonlinear_find_concentrations_witherror(inten_test_sliced, A_sliced, B_sliced, ...
            @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
            alpha, beta, maxIterCnt, eps, 10^alpha_C, 1));

    A_all{i} = A;
    A_sliced_all{i} = A_sliced;
    B_all{i} = B;
    B_sliced_all{i} = B_sliced;
    C_train_all{i} = C_train;
    C_test_all{i} = C_test;
    
    save('sample_size_frankenshtein_fixed_reg_m5', 'train_error', 'test_error', ...
        'A_all', 'A_sliced_all', 'B_all', 'B_sliced_all', 'C_train_all', ...
        'C_test_all');
    
%     [A_train, C_train, ~, A_train_sliced] = calibrate_model_parallel(inten_train, ...
%         inten_train_sliced, inten_full_idx_curr, @(I) nmf_rma(I));
%     [A_test_1, C_test_1, ~, A_test_sliced_1] = calibrate_model_parallel(inten_test, ...
%         inten_test_sliced, inten_full_idx_curr, @(I) nmf_rma(I));
%     [C_test_2] = find_concentrations(inten_test_sliced, A_train_sliced, @(I, A) nmf_rma_fixedA(I, A));
% 
%     [A_train_LO, C_train_LO, ~, A_train_sliced_LO] = calibrate_model_parallel(inten_train, inten_train_sliced, ...
%         inten_full_idx_curr, @(I) nmf_rma_LO_reg(I));
%     [C_test_LO] = find_concentrations(inten_test_sliced, ...
%         A_train_sliced_LO, @(I, A) nmf_rma_LO_reg_fixedA(I, A));
% 
%     train_error_rma(i) = quality_functional_linear(inten_train_sliced, A_train_sliced, C_train);
%     test_error_rma_1(i) = quality_functional_linear(inten_test_sliced, A_test_sliced_1, C_test_1);
%     test_error_rma_2(i) = quality_functional_linear(inten_test_sliced, A_train_sliced, C_test_2);
%     test_error_rma_LO(i) = quality_functional_linear(inten_test_sliced, A_train_sliced_LO, C_test_LO);
%     
%     save('sample_size_rma_1', 'train_error_rma', 'test_error_rma_1', 'test_error_rma_2', ...
%         'test_error_rma_LO');
end
        