alpha = 2;
beta = 1;
maxIterCnt = 1000;
eps = 1e-5;

alpha_C_range = -50:1:20;

train_error = zeros(length(alpha_C_range), 1);
validation_error = zeros(length(alpha_C_range), 1);
test_error = zeros(length(alpha_C_range), 1);

A_all = cell(length(alpha_C_range), 1);
A_sliced_all = cell(length(alpha_C_range), 1);
B_all = cell(length(alpha_C_range), 1);
B_sliced_all = cell(length(alpha_C_range), 1);
C_validation_all = cell(length(alpha_C_range), 1);
C_train_all = cell(length(alpha_C_range), 1);
C_test_all = cell(length(alpha_C_range), 1);

arrays_cnt = size(inten_full, 2) - 2;
%validation_arrays = randsample(arrays_cnt, 200);
%test_arrays = randsample(setdiff(1:arrays_cnt, validation_arrays), 500);
load('test_arrays_idx');

probes_cnt = inten_full_idx{50}(end);
inten_validation = inten_full(1:probes_cnt, [validation_arrays' end-1 end]);
inten_test = inten_full(1:probes_cnt, [test_arrays end-1 end]);
inten_validation_sliced = inten_full_sliced(1:50);
inten_test_sliced = inten_full_sliced(1:50);
for i = 1:50
    inten_validation_sliced{i} = inten_validation_sliced{i}(:, validation_arrays);
    inten_test_sliced{i} = inten_test_sliced{i}(:, test_arrays);
end

inten_full_idx_curr = inten_full_idx(1:50);

alpha_C = 1;

inten_full_curr = inten_full(1:probes_cnt, [setdiff(1:arrays_cnt, union(validation_arrays, test_arrays)) end-1 end]);
inten_full_sliced_curr = inten_full_sliced(1:50);
for i = 1:50
    inten_full_sliced_curr{i} = inten_full_sliced_curr{i}(:, setdiff(1:arrays_cnt, union(validation_arrays, test_arrays)));
end

train_size = 300;

inten_train = zeros(probes_cnt, train_size+2);
inten_train_sliced = inten_full_sliced_curr;
arrays_train = cell(50, 1);
for j = 1:50
    arrays_train{j} = pick_arrays(inten_full_sliced_curr{j}, train_size);
    inten_train_sliced{j} = inten_full_sliced_curr{j}(:, arrays_train{j});
    inten_train(inten_full_idx{j}, :) = inten_full_curr(inten_full_idx{j}, [arrays_train{j} end-1 end]);
end

for i = 1:length(alpha_C_range)
    alpha_C = alpha_C_range(i);
    
    fprintf('%d\n', alpha_C);
    
    [A, B, C_train, ~, ~, A_sliced, B_sliced] = nonlinear_calibrate_model_short(inten_train, ...
        inten_train_sliced, inten_full_idx_curr, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_C, 1));
    
    [C_validation, test_error(i)] = nonlinear_find_concentrations_witherror(inten_validation_sliced, A_sliced, B_sliced, ...
            @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
            alpha, beta, maxIterCnt, eps, 10^alpha_C, 1));
    
    [C_test, test_error(i)] = nonlinear_find_concentrations_witherror(inten_test_sliced, A_sliced, B_sliced, ...
            @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
            alpha, beta, maxIterCnt, eps, 10^alpha_C, 1));
        
    train_error(i) = quality_functional(inten_train_sliced, A_sliced, B_sliced, C_train);
    validation_error(i) = quality_functional(inten_validation_sliced, A_sliced, B_sliced, C_validation);
    test_error(i) = quality_functional(inten_test_sliced, A_sliced, B_sliced, C_test);
    
    A_all{i} = A;
    A_sliced_all{i} = A_sliced;
    B_all{i} = B;
    B_sliced_all{i} = B_sliced;
    C_train_all{i} = C_train;
    C_validation_all{i} = C_validation;
    C_test_all{i} = C_test;
    
    save('sample_size_frankenshtein_all_reg', 'train_error', 'test_error', 'validation_error', ...
        'A_all', 'A_sliced_all', 'B_all', 'B_sliced_all', 'C_train_all', ...
        'C_validation_all', 'C_test_all', 'arrays_train', 'validation_arrays', ...
        'test_arrays');
end
        