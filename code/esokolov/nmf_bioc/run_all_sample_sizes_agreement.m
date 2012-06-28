full_sample_size = 1000;
sample_size_range = 50:50:1000;
factorize_smart = @(I) nmf_smart(I, 1, 1000, 1e-6);

quality_arrays_mad_smart = zeros(length(sample_size_range), 1);
quality_arrays_ouliers_smart = zeros(length(sample_size_range), 1);
quality_probes_mad_smart = zeros(length(sample_size_range), 1);
quality_probes_ouliers_smart = zeros(length(sample_size_range), 1);

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
    [partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl] = ...
        generate_partitions(inten, inten_sliced, inten_full_idx);
    fprintf(' Done\n');
    
    fprintf('Factorizing I...');
    [quality_arrays_mad_smart(i), quality_arrays_ouliers_smart(i), quality_probes_mad_smart(i), quality_probes_ouliers_smart(i), ...
        arrays_factors_smart, probes_factors_smart] = ...
        get_quality_cv(inten, inten_sliced, inten_full_idx, ...
        factorize_smart, partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl);
    fprintf(' Done\n');
    
    fprintf('Mad(A1 - A2): %f\nMad(C1 - C2): %f\n\n', quality_arrays_mad_smart(i), quality_probes_mad_smart(i));
    
    save(['sample_size_agr_1000_' int2str(sample_size)], 'quality_arrays_mad_smart', ...
        'quality_arrays_ouliers_smart', 'quality_probes_mad_smart', 'quality_probes_ouliers_smart', ...
        'arrays_factors_smart', 'probes_factors_smart');
end