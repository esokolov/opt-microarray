full_sample_size = 900;
sample_size_range = 50:50:900;

inten_test = inten_full(:, 901:end);
inten_test_sliced = inten_full_sliced;
for i = 1:length(inten_test_sliced)
    inten_test_sliced{i} = inten_test_sliced{i}(:, 901:end);
end

train_error = zeros(length(sample_size_range), 1);
test_error = zeros(length(sample_size_range), 1);

alpha = -0.5;
beta = 1;

for i = 1:length(sample_size_range)
    sample_size = sample_size_range(i);
    fprintf('Training sample size: %d\n', sample_size);
    
    subsample = randsample(full_sample_size, sample_size, false)';
    inten = inten_full(:, [subsample size(inten_full, 2)-1:size(inten_full, 2)]);
    inten_sliced = inten_full_sliced;
    for j = 1:length(inten_sliced)
        inten_sliced{j} = inten_sliced{j}(:, subsample);
    end
    
    fprintf('Factorizing I...');
    [A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, inten_full_idx, ...
        @(I) nmf_alpha_beta(I, 1, alpha, beta, 50, 1e-6));
    fprintf(' Done\n');
    
    %C_test = max(eps, (A' * A) \ (A' * inten_test(:, 1:100)));
    fprintf('Finding C for testing sample...');
    C_test = find_concentrations(inten_test_sliced, A_sliced, ...
        @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, 50, 1e-6));
    fprintf(' Done\n');
    
	train_error(i) = nmf_alpha_beta_divergence(inten(:, 1:sample_size), A * C, alpha, beta);
    test_error(i) = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta);
    
    fprintf('Train error: %e\nTest error: %e\n\n', train_error(i), test_error(i));
    
    save(['sample_size_900_' int2str(sample_size)], 'train_error', 'test_error');
end