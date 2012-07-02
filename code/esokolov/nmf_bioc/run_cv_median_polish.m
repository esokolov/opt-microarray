
inten = inten_full(:, [1:100 1001:1002]);
inten_sliced = inten_full_sliced;
for i = 1:length(inten_sliced)
    inten_sliced{i} = inten_sliced{i}(:, 1:100);
end

inten_test = inten_full(:, 901:end);
inten_test_sliced = inten_full_sliced;
for i = 1:length(inten_test_sliced)
    inten_test_sliced{i} = inten_test_sliced{i}(:, 901:end);
end

%fprintf('Generating partitions...');
%[partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
%   generate_partitions(inten, inten_sliced, inten_full_idx);
%fprintf(' Done\n');
load('alpha_beta_partition.mat');

fprintf('Factorizing partitions of I...');
[quality_arrays_mad_mp, quality_arrays_ouliers_mp, quality_probes_mad_mp, quality_probes_ouliers_mp, ...
    arrays_factors_smart, probes_factors_smart, converged_cnt_arrays1_mp, converged_cnt_arrays2_mp, ...
    converged_cnt_probes1_mp, converged_cnt_probes2_mp] = ...
    get_quality_cv(inten, inten_sliced, inten_full_idx, ...
    @(I) nmf_median_polish(I, 50, 1e-6), partition_arrays, partition_probes, partition_probes_compl, ...
    partition_probes_sliced, partition_probes_sliced_compl);
fprintf(' Done\n');

conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
for k = 1:length(conc_gene_dist)
    conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
end
concentration_dists_mp = mad(conc_gene_dist);

goodness_of_fit_mp = norm(inten(:, partition_arrays) - arrays_factors_smart{1} * arrays_factors_smart{2}, 'fro');

% fprintf('Factorizing I...');
% [A, C, Avect, A_sliced, converged_cnt_full_mp] = calibrate_model_parallel(inten, inten_sliced, ...
%     inten_full_idx, @(I) nmf_median_polish(I, 50, 1e-6));
% fprintf(' Done\n');
% 
% fprintf('Factorizing I_test with fixed A...');
% [C_test converged_cnt_fixedA_mp] = find_concentrations(inten_test_sliced, A_sliced, ...
%     @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, 50, 1e-6));
% fprintf(' Done\n');
% 
% overfitting_native_mp = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta) - ...
%     nmf_alpha_beta_divergence(inten(:, 1:end-2), A * C, alpha, beta);
% overfitting_fro_mp = norm(inten_test(:, 1:end-2) - A * C_test, 'fro') - ...
%     norm(inten(:, 1:end-2) - A * C, 'fro');

fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %e\nMad(dist(C1, C2)): %e\nGoodness of fit: %e\n\n', ...
    quality_arrays_mad_mp, quality_probes_mad_mp, concentration_dists_mp, goodness_of_fit_mp);

save('alpha_beta_cv_median_polish', 'goodness_of_fit_mp', 'quality_arrays_mad_mp', ...
    'concentration_dists_mp', ...
    'quality_arrays_ouliers_mp', 'quality_probes_mad_mp', 'quality_probes_ouliers_mp', ...
    'converged_cnt_arrays1_mp', 'converged_cnt_arrays2_mp', 'converged_cnt_probes1_mp', 'converged_cnt_probes2_mp');%, ...
    %'arrays_factors_smart', 'probes_factors_smart');