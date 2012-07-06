
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
[quality_arrays_mad_rma, quality_arrays_ouliers_rma, quality_probes_mad_rma, quality_probes_ouliers_rma, ...
    arrays_factors_smart, probes_factors_smart, converged_cnt_arrays1_rma, converged_cnt_arrays2_rma, ...
    converged_cnt_probes1_rma, converged_cnt_probes2_rma] = ...
    get_quality_cv(inten, inten_sliced, inten_full_idx, ...
    @(I) nmf_rma(I), partition_arrays, partition_probes, partition_probes_compl, ...
    partition_probes_sliced, partition_probes_sliced_compl);
fprintf(' Done\n');

conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
for k = 1:length(conc_gene_dist)
    conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
end
concentration_dists_rma = mad(conc_gene_dist);

goodness_of_fit_rma = norm(inten(:, partition_arrays) - arrays_factors_smart{1} * arrays_factors_smart{2}, 'fro');

% fprintf('Factorizing I...');
% [A, C, Avect, A_sliced, converged_cnt_full_rma] = calibrate_model_parallel(inten, inten_sliced, ...
%     inten_full_idx, @(I) nmf_median_polish(I, 1000, 1e-6));
% fprintf(' Done\n');
% 
% fprintf('Factorizing I_test with fixed A...');
% [C_test converged_cnt_fixedA_rma] = find_concentrations(inten_test_sliced, A_sliced, ...
%     @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, 50, 1e-6));
% fprintf(' Done\n');
% 
% overfitting_native_rma = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta) - ...
%     nmf_alpha_beta_divergence(inten(:, 1:end-2), A * C, alpha, beta);
% overfitting_fro_rma = norm(inten_test(:, 1:end-2) - A * C_test, 'fro') - ...
%     norm(inten(:, 1:end-2) - A * C, 'fro');

fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %e\nMad(dist(C1, C2)): %e\nGoodness of fit: %e\n\n', ...
    quality_arrays_mad_rma, quality_probes_mad_rma, concentration_dists_rma, goodness_of_fit_rma);

save('alpha_beta_cv_median_polish_supernew', 'goodness_of_fit_rma', 'quality_arrays_mad_rma', ...
    'concentration_dists_rma', ...
    'quality_arrays_ouliers_rma', 'quality_probes_mad_rma', 'quality_probes_ouliers_rma', ...
    'converged_cnt_arrays1_rma', 'converged_cnt_arrays2_rma', 'converged_cnt_probes1_rma', 'converged_cnt_probes2_rma');%, ...
    %'arrays_factors_smart', 'probes_factors_smart');