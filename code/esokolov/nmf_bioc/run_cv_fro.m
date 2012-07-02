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

fprintf('Factorizing partitions of I...');
[quality_arrays_mad_fro, quality_arrays_ouliers_fro, quality_probes_mad_fro, quality_probes_ouliers_fro, ...
    arrays_factors_smart, probes_factors_smart] = ...
    get_quality_cv(inten, inten_sliced, inten_full_idx, ...
    @(I) nnmf(I, 1), partition_arrays, partition_probes, partition_probes_compl, ...
    partition_probes_sliced, partition_probes_sliced_compl);
fprintf(' Done\n');

conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
for k = 1:length(conc_gene_dist)
    conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
end
concentration_dists_fro = mad(conc_gene_dist);

goodness_of_fit_fro = norm(inten(:, partition_arrays) - arrays_factors_smart{1} * arrays_factors_smart{2}, 'fro');

fprintf('Factorizing I...');
[A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, ...
    inten_full_idx, @(I) nnmf(I, 1));
fprintf(' Done\n');

fprintf('Factorizing I_test with fixed A...');
C_test = find_concentrations(inten_test_sliced, A_sliced, ...
    @(I_arg, A_arg) max(eps, pinv(min(full(A_arg' * A_arg), 1e150)) * (A_arg' * I_arg)));
%C_test = max(eps, pinv(min(full(A' * A), 1e150)) * (A' * I));
fprintf(' Done\n');

overfitting_native_fro = norm(inten_test(:, 1:end-2) - A * C_test, 'fro') - ...
    norm(inten(:, 1:end-2) - A * C, 'fro');
overfitting_fro_fro = norm(inten_test(:, 1:end-2) - A * C_test, 'fro') - ...
    norm(inten(:, 1:end-2) - A * C, 'fro');

fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %e\nMad(dist(C1, C2)): %e\nGoodness of fit: %e\nOverfitting (alpha-beta): %e\nOverfitting (L2): %e\n\n', ...
    quality_arrays_mad_fro, quality_probes_mad_fro, concentration_dists_fro, goodness_of_fit_fro, ...
    overfitting_native_fro, overfitting_fro_fro);

save('alpha_beta_cv_fro.mat', 'goodness_of_fit_fro', 'quality_arrays_mad_fro', ...
    'overfitting_native_fro', 'overfitting_fro_fro', 'concentration_dists_fro', ...
    'quality_arrays_ouliers_fro', 'quality_probes_mad_fro', 'quality_probes_ouliers_fro');%, ...
    %'arrays_factors_smart', 'probes_factors_smart');