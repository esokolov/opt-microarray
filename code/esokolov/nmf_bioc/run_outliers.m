training_sample_size = 2500;
last_idx = size(inten_full, 2);

inten = inten_full(:, [1:training_sample_size last_idx-1:last_idx]);
inten_sliced = inten_full_sliced;
for i = 1:length(inten_sliced)
    inten_sliced{i} = inten_sliced{i}(:, 1:training_sample_size);
end

inten_test = inten_full(:, training_sample_size+1:end);
inten_test_sliced = inten_full_sliced;
for i = 1:length(inten_test_sliced)
    inten_test_sliced{i} = inten_test_sliced{i}(:, training_sample_size+1:end);
end

fprintf('Generating partitions...');
[partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
   generate_partitions(inten, inten_sliced, inten_full_idx);
fprintf(' Done\n');

alpha = -0.5;
beta = 1;
maxIterCnt = 500;



[quality_arrays_mad_full, quality_arrays_ouliers_full, quality_probes_mad_full, quality_probes_ouliers_full, ...
    arrays_factors_smart, probes_factors_smart] = ...
    get_quality_cv(inten, inten_sliced, inten_full_idx, ...
    @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6), partition_arrays, partition_probes, partition_probes_compl, ...
    partition_probes_sliced, partition_probes_sliced_compl);

conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
for k = 1:length(conc_gene_dist)
    conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
end
concentration_dist_full = mad(conc_gene_dist);

fprintf('Arrays quality full: %e\nProbes quality full: %e\nConcentration dist full: %e\n\n', ...
    quality_arrays_mad_full, quality_probes_mad_full, concentration_dist_full);

filtered_cnt = 100:100:training_sample_size;
quality_arrays_filtered = zeros(length(filtered_cnt), 1);
quality_probes_filtered = zeros(length(filtered_cnt), 1);
concentration_dists_filtered = zeros(length(filtered_cnt), 1);

for i = 1:length(filtered_cnt)
    fprintf('%d\n', i);
    
    inten_curr = inten(:, [1:filtered_cnt(i) size(inten, 2)-1:size(inten, 2)]);
    inten_curr_sliced = inten_sliced;
    for j = 1:length(inten_curr_sliced)
        inten_curr_sliced{j} = inten_curr_sliced{j}(:, 1:filtered_cnt(i));
    end
    
    partition_arrays = generate_partitions(inten_curr, inten_curr_sliced, inten_full_idx);
    
    [quality_arrays_filtered(i), ~, quality_probes_filtered(i), ~, ...
        arrays_factors_smart, probes_factors_smart] = ...
        get_quality_cv(inten_curr, inten_curr_sliced, inten_full_idx, ...
        @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6), partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl);
        
    conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
    for k = 1:length(conc_gene_dist)
        conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
    end
    concentration_dists_filtered(i) = mad(conc_gene_dist);
    
    fprintf('Arrays quality: %e\nProbes quality: %e\nConcentration dist: %e\n\n', ...
    quality_arrays_filtered(i), quality_probes_filtered(i), concentration_dists_filtered(i));
end




[A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, ...
            inten_full_idx, @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6));
        
C_test = find_concentrations(inten_test_sliced, A_sliced, ...
            @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, maxIterCnt, 1e-6));

err_full = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta);
fprintf('Quality with full training sample: %e\n', err_full);
        
AC = A * C;
dists = zeros(training_sample_size, 1);
for i = 1:length(dists)
    dists(i) = nmf_alpha_beta_divergence(inten(:, i), AC(:, i), alpha, beta);
end

[~, idx] = sort(dists);
dists = dists(idx);
inten = inten(:, [idx' size(inten, 2)-1:size(inten, 2)]);

filtered_cnt = 100:100:training_sample_size;
err_filtered = zeros(length(filtered_cnt), 1);

for i = 1:length(filtered_cnt)
    fprintf('%d\n', i);
    
    inten_curr = inten(:, [1:filtered_cnt(i) size(inten, 2)-1:size(inten, 2)]);
    inten_curr_sliced = inten_sliced;
    for j = 1:length(inten_curr_sliced)
        inten_curr_sliced{j} = inten_curr_sliced{j}(:, 1:filtered_cnt(i));
    end
    
    [A, C, Avect, A_sliced] = calibrate_model_parallel(inten_curr, inten_curr_sliced, ...
            inten_full_idx, @(I) nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, 1e-6));
        
    C_test = find_concentrations(inten_test_sliced, A_sliced, ...
                @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, maxIterCnt, 1e-6));

    err_filtered(i) = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta);
    
    fprintf('Quality with %d samples: %e\n', filtered_cnt(i), err_filtered(i));
end

save('outliers_res.mat', 'err_full', 'err_filtered', ...
    'quality_arrays_mad_full', 'quality_probes_mad_full', 'concentration_dist_full', ...
    'quality_arrays_filtered', 'quality_probes_filtered', 'concentration_dists_filtered');
