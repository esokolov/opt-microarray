full_sample_size = 1000;
sample_size_range = 50:50:1000;
%factorize_smart = @(I) nmf_smart(I, 1, 1000, 1e-6);

quality_arrays_mad = zeros(length(sample_size_range), 1);
quality_arrays_ouliers = zeros(length(sample_size_range), 1);
quality_probes_mad = zeros(length(sample_size_range), 1);
quality_probes_ouliers = zeros(length(sample_size_range), 1);
concentration_dists = zeros(length(sample_size_range), 1);

load('alpha_beta_partition.mat');

for i = 1:length(sample_size_range)
    sample_size = sample_size_range(i);
    fprintf('Training sample size: %d\n', sample_size);
    
    subsample = randsample(full_sample_size, sample_size, false)';
    inten = inten_full(:, [subsample size(inten_full, 2)-1:size(inten_full, 2)]);
    inten_sliced = inten_full_sliced;
    for j = 1:length(inten_sliced)
        inten_sliced{j} = inten_sliced{j}(:, subsample);
    end
    
    fprintf('Generating partitions...');
    [partition_arrays, ~, ~, ~, ~] = ...
        generate_partitions(inten, inten_sliced, inten_full_idx);
    fprintf(' Done\n');
    
    fprintf('Factorizing I...');
    [quality_arrays_mad(i), quality_arrays_ouliers(i), quality_probes_mad(i), quality_probes_ouliers(i), ...
        arrays_factors_smart, probes_factors_smart] = ...
        get_quality_cv(inten, inten_sliced, inten_full_idx, ...
        @(I) nmf_alpha_beta(I, 1, -1, 1, 50, 1e-6), partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl);
    fprintf(' Done\n');
    
    conc_gene_dist = zeros(size(probes_factors_smart{2}, 1), 1);
    for k = 1:length(conc_gene_dist)
        conc_gene_dist(k) = Cdist(probes_factors_smart{2}(k, :), probes_factors_smart{6}(k, :));
    end
    concentration_dists(i) = mad(conc_gene_dist);
    
    fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %f\nc_dist: %f\n\n', quality_arrays_mad(i), quality_probes_mad(i), concentration_dists(i));
    
    save(['sample_size_agr_900_' int2str(sample_size)], 'quality_arrays_mad', ...
        'quality_arrays_ouliers', 'quality_probes_mad', 'quality_probes_ouliers', 'concentration_dists');
end