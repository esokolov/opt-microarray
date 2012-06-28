full_sample_size = 500;
sample_size_range = [10:25:500 500];

train_error = zeros(length(sample_size_range), 1);
test_error = zeros(length(sample_size_range), 1);


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
    [A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, inten_full_idx, factorize_smart);
    fprintf(' Done\n');
    
    %C_test = max(eps, (A' * A) \ (A' * inten_test(:, 1:100)));
    fprintf('Finding C for testing sample...');
    C_test = find_concentrations(inten_test_sliced, A_sliced, factorize_smart_fixedA);
    fprintf(' Done\n');
    
	train_error(i) = nmf_alpha_divergence(inten(:, 1:sample_size), A * C, 0);
    test_error(i) = nmf_alpha_divergence(inten_test(:, 1:end-2), A * C_test, 0);
    
    fprintf('Train error: %e\nTest error: %e\n\n', train_error(i), test_error(i));
    
    save(['sample_size_500_' int2str(sample_size)], 'A', 'C', 'Avect', 'A_sliced', 'C_test', 'train_error', 'test_error');
end