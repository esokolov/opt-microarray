alpha_range = -3:0.5:3;
beta_range = -2:0.5:10;

quality_arrays_mad = zeros(length(alpha_range), length(beta_range));
quality_arrays_ouliers = zeros(length(alpha_range), length(beta_range));
quality_probes_mad = zeros(length(alpha_range), length(beta_range));
quality_probes_ouliers = zeros(length(alpha_range), length(beta_range));
concentration_dists = zeros(length(alpha_range), length(beta_range));
goodness_of_fit = zeros(length(alpha_range), length(beta_range));
overfitting_native = zeros(length(alpha_range), length(beta_range));
overfitting_fro = zeros(length(alpha_range), length(beta_range));

%inten = inten_full;
%inten_sliced = inten_full_sliced;

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

fprintf('Generating partitions...');
[partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
   generate_partitions(inten, inten_sliced, inten_full_idx);
fprintf(' Done\n');

for i = 1:length(alpha_range)
    for j = 1:length(beta_range)
        alpha = alpha_range(i);
        beta = beta_range(j);
        fprintf('Alpha = %f, Beta = %f\n', alpha, beta);

        fprintf('Factorizing partitions of I...');
        [quality_arrays_mad(i, j), quality_arrays_ouliers(i, j), quality_probes_mad(i, j), quality_probes_ouliers(i, j), ...
            arrays_factors_smart, probes_factors_smart] = ...
            get_quality_cv(inten, inten_sliced, inten_full_idx, ...
            @(I) nmf_alpha_beta(I, 1, alpha, beta, 50, 1e-6), partition_arrays, partition_probes, partition_probes_compl, ...
            partition_probes_sliced, partition_probes_sliced_compl);
        fprintf(' Done\n');
        
        conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
        for k = 1:length(conc_gene_dist)
            conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
        end
        concentration_dists(i, j) = mad(conc_gene_dist);
        
        goodness_of_fit(i, j) = norm(inten(:, partition_arrays) - arrays_factors_smart{1} * arrays_factors_smart{2}, 'fro');
        
        fprintf('Factorizing I...');
        [A, C, Avect, A_sliced] = calibrate_model_parallel(inten, inten_sliced, ...
            inten_full_idx, @(I) nmf_alpha_beta(I, 1, alpha, beta, 50, 1e-6));
        fprintf(' Done\n');
        
        fprintf('Factorizing I_test with fixed A...');
        C_test = find_concentrations(inten_test_sliced, A_sliced, ...
            @(I_arg, A_arg) nmf_alpha_beta_fixedA(I_arg, A_arg, alpha, beta, 50, 1e-6));
        fprintf(' Done\n');
        
        overfitting_native(i, j) = nmf_alpha_beta_divergence(inten_test(:, 1:end-2), A * C_test, alpha, beta) - ...
            nmf_alpha_beta_divergence(inten(:, 1:end-2), A * C, alpha, beta);
        overfitting_fro(i, j) = norm(inten_test(:, 1:end-2) - A * C_test, 'fro') - ...
            norm(inten(:, 1:end-2) - A * C, 'fro');

        fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %e\nMad(dist(C1, C2)): %e\nGoodness of fit: %e\nOverfitting (alpha-beta): %e\nOverfitting (L2): %e\n\n', ...
            quality_arrays_mad(i, j), quality_probes_mad(i, j), concentration_dists(i, j), goodness_of_fit(i, j), ...
            overfitting_native(i, j), overfitting_fro(i, j));

        save(['alpha_beta_cv_second_' num2str(alpha - min(alpha_range)) '_' num2str(beta - min(beta_range))], 'goodness_of_fit', 'quality_arrays_mad', ...
            'overfitting_native', 'overfitting_fro', 'concentration_dists', ...
            'quality_arrays_ouliers', 'quality_probes_mad', 'quality_probes_ouliers');%, ...
            %'arrays_factors_smart', 'probes_factors_smart');
    end
end