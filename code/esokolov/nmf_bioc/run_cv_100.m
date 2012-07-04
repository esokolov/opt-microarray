iterCnt = 20;

full_sample_size = 1000;

train_error = zeros(iterCnt, 1);
test_error = zeros(iterCnt, 1);

alpha = 1;
beta = -0.5;

for i = 1:iterCnt
    fprintf('Iteration %d\n', i);
    
    subsample = randsample(full_sample_size, 200, false)';
    
    inten = inten_full(:, [subsample(1:100) size(inten_full, 2)-1:size(inten_full, 2)]);
    inten_sliced = inten_full_sliced;
    for j = 1:length(inten_sliced)
        inten_sliced{j} = inten_sliced{j}(:, subsample(1:100));
    end
    
    inten_test = inten_full(:, [subsample(101:200) size(inten_full, 2)-1:size(inten_full, 2)]);
    inten_test_sliced = inten_full_sliced;
    for j = 1:length(inten_sliced)
        inten_test_sliced{j} = inten_test_sliced{j}(:, subsample(101:200));
    end
    
    fprintf('Factorizing I...');
    [A, C, Avect, A_sliced, cc] = calibrate_model_parallel(inten, inten_sliced, inten_full_idx, ...
        @(I) nmf_alpha_beta(I, 1, alpha, beta, 50, 1e-6));
    fprintf(' Done\n');
    
    %C_test = max(eps, (A' * A) \ (A' * inten_test(:, 1:100)));
    fprintf('Finding C for testing sample...');
    C_test = find_concentrations(inten_test_sliced, A_sliced, ...
        @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, 50, 1e-6));
    fprintf(' Done\n');
    
	train_error(i) = nmf_alpha_beta_divergence(inten(:, 1:end-2), A * C, alpha, beta);
    test_error(i) = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta);
    
    fprintf('Train error: %e\nTest error: %e\n\n', train_error(i), test_error(i));
    
    if (test_error(i) < train_error(i))
        break;
    end
    
    save('cv_100.mat', 'train_error', 'test_error');
end