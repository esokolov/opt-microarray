train_sizes = 100:100:1000;
alpha = 2;
beta = 1;
maxIterCnt = 1000;
eps = 1e-5;

reg_best = zeros(length(train_sizes), 1);
train_error = zeros(length(train_sizes), 1);
validation_error = zeros(length(train_sizes), 1);
test_error = zeros(length(train_sizes), 1);

A_all = cell(length(train_sizes), 1);
A_sliced_all = cell(length(train_sizes), 1);
B_all = cell(length(train_sizes), 1);
B_sliced_all = cell(length(train_sizes), 1);
C_train_all = cell(length(train_sizes), 1);
C_validation_all = cell(length(train_sizes), 1);
C_test_all = cell(length(train_sizes), 1);

for i = 1:length(train_sizes)
    train_size = train_sizes(i)
    
    inten_train_tmp = inten_train(:,[1:train_size end-1 end]);
    inten_train_sliced_tmp = inten_train_sliced;
    for j = 1:50
        inten_train_sliced_tmp{j} = inten_train_sliced{j}(:, 1:train_size);
    end

    [A, B, C_train, C_validation, A_sliced, B_sliced, reg_best(i), train_error(i), validation_error(i)] = ...
            nonlinear_calibrate_frankenstein_model_exh_witherror(inten_train_tmp, inten_train_sliced_tmp, inten_validation_sliced, inten_full_idx, ...
            alpha, beta, maxIterCnt, eps);
        
%    train_error(i) = quality_functional(inten_train_sliced_tmp, A_sliced, B_sliced, C_train);
%    validation_error(i) = quality_functional(inten_validation_sliced, A_sliced, B_sliced, C_validation);
    
    [C_test, test_error(i)] = nonlinear_find_concentrations_witherror(inten_test_sliced, A_sliced, B_sliced, ...
        @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
        alpha, beta, maxIterCnt, eps, reg_best(i), 1));
    
%    test_error(i) = quality_functional(inten_test_sliced, A_sliced, B_sliced, C_test);
    
    A_all{i} = A;
    A_sliced_all{i} = A_sliced;
    B_all{i} = B;
    B_sliced_all{i} = B_sliced;
    C_train_all{i} = C_train;
    C_validation_all{i} = C_validation;
    C_test_all{i} = C_test;
    
    save('sample_size_frankenshtein_exh.mat', 'train_error', 'test_error', 'validation_error', ...
        'A_all', 'A_sliced_all', 'B_all', 'B_sliced_all', 'C_train_all', ...
        'C_validation_all', 'C_test_all', 'reg_best', 'arrays_train', 'validation_arrays', ...
        'test_arrays');
end        