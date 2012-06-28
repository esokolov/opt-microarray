%% функции ошибки
factorize_frobenius = @(I) nnmf(I, 1);
factorize_smart = @(I) nmf_smart(I, 1, 1000, 1e-6);
factorize_KL = @(I) nmf_alpha(I, 1, 1, 1000, 1e-6);
factorize_hellinger = @(I) nmf_alpha(I, 1, 0.5, 1000, 1e-6);
factorize_chi_squared = @(I) nmf_alpha(I, 1, 2, 1000, 1e-6);
factorize_be = @(I) nmf_smart_generalized(I, 1, 1000, 1e-6, 'BE1');
factorize_rjs = @(I) nmf_smart_generalized(I, 1, 1000, 1e-6, 'RJS');

factorize_smart_fixedA = @(I, A) nmf_smart_fixedA(I, A, 1000, 1e-6);

%% разбиения
sampleSize = size(inten, 2) - 2;
partitionSize = fix(sampleSize / 2);
%partition_arrays = datasample(1:sampleSize, partitionSize, 'Replace', false);
partition_arrays = randsample(1:sampleSize, partitionSize, false);

probesets = unique(inten(:, end));
partition_probes = zeros(size(inten, 1), 1);
partition_probes_compl = zeros(size(inten, 1), 1);
partitionSize = 0;
partitionComplSize = 0;
for i = 1:length(probesets)
    idx = find(inten(:, end) == probesets(i));
    probesetSize = length(idx);
    if (probesetSize == 1)
        continue;
    end
    currPartitionSize = fix(probesetSize / 2);
    
    currPartition = randsample(idx, currPartitionSize, false);
    partition_probes((partitionSize + 1):(partitionSize + currPartitionSize)) = currPartition;
    partitionSize = partitionSize + currPartitionSize;
    
    currPartitionCompl = setdiff(idx, currPartition);
    partition_probes_compl((partitionComplSize + 1):(partitionComplSize + length(currPartitionCompl))) = ...
        currPartitionCompl;
    partitionComplSize = partitionComplSize + length(currPartitionCompl);
end
partition_probes = partition_probes(1:partitionSize);
partition_probes_compl = partition_probes_compl(1:partitionComplSize);

%%
[quality_arrays_mad_fro, quality_arrays_ouliers_fro, quality_probes_mad_fro, quality_probes_ouliers_fro, ...
        arrays_factors_fro, probes_factors_fro] = ...
        get_quality_cv(inten, factorize_frobenius, partition_arrays, partition_probes, partition_probes_compl);
[quality_arrays_mad_smart, quality_arrays_ouliers_smart, quality_probes_mad_smart, quality_probes_ouliers_smart, ...
        arrays_factors_smart, probes_factors_smart] = ...
        get_quality_cv(inten, factorize_smart, partition_arrays, partition_probes, partition_probes_compl);
[quality_arrays_mad_KL, quality_arrays_ouliers_KL, quality_probes_mad_KL, quality_probes_ouliers_KL, ...
        arrays_factors_KL, probes_factors_KL] = ...
        get_quality_cv(inten, factorize_KL, partition_arrays, partition_probes, partition_probes_compl);
[quality_arrays_mad_hellinger, quality_arrays_ouliers_hellinger, quality_probes_mad_hellinger, quality_probes_ouliers_hellinger, ...
        arrays_factors_hellinger, probes_factors_hellinger] = ...
        get_quality_cv(inten, factorize_hellinger, partition_arrays, partition_probes, partition_probes_compl);
[quality_arrays_mad_chi_squared, quality_arrays_ouliers_chi_squared, quality_probes_mad_chi_squared, quality_probes_ouliers_chi_squared, ...
        arrays_factors_chi_squared, probes_factors_chi_squared] = ...
        get_quality_cv(inten, factorize_chi_squared, partition_arrays, partition_probes, partition_probes_compl);
%%
[quality_arrays_mad_be, quality_arrays_ouliers_be, quality_probes_mad_be, quality_probes_ouliers_be, ...
        arrays_factors_be, probes_factors_be] = ...
        get_quality_cv(inten, factorize_be, partition_arrays, partition_probes, partition_probes_compl);
[quality_arrays_mad_rjs, quality_arrays_ouliers_rjs, quality_probes_mad_rjs, quality_probes_ouliers_rjs, ...
        arrays_factors_rjs, probes_factors_rjs] = ...
        get_quality_cv(inten, factorize_rjs, partition_arrays, partition_probes, partition_probes_compl);
    
%% графики качества на отдельных чипах
[A1, C1] = calibrate_model(inten, factorize_KL);
[A2, C2] = calibrate_model(inten, factorize_hellinger);
%%
I1 = A1 * C1;
I2 = A2 * C2;
I1(isnan(I1)) = 0;
I2(isnan(I2)) = 0;
q1 = zeros(51, 1);
q2 = zeros(51, 1);
for i = 1:51
    q1(i) = nmf_alpha_divergence(inten(:, i), I1(:, i), 1);
    q2(i) = nmf_alpha_divergence(inten(:, i), I2(:, i), 0.5);
end
%%
plot(q1, q2, 'o');

%%
[A_smart_full, C_smart_full, Avect_smart_full] = calibrate_model(intens, factorize_smart);

%% одна оптимизация, все пары критериев
cost_functions = { (@(I, AC) norm(I - AC, 'fro'));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 1));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0.5));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 2)); };
costs_all = zeros(100, 5);
I_hat = A_smart_full * C_smart_full;
for i = 1:5
    for j = 1:100
        costs_all(j, i) = cost_functions{i}(intens(:, j), I_hat(:, j));
    end
end

%% все расстояния между концентрациями
cost_functions = { (@(I, AC) norm(I - AC, 'fro'));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 1));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0.5));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 2));
                   (@(I, AC) nmf_alpha_divergence(I, AC, -1)); };
costs_all = zeros(100, 6);
%I_hat = A_smart_full * C_smart_full;
for i = 1:6
    for j = 1:100
        costs_all(j, i) = cost_functions{i}(probes_factors_smart{2}(:, j), probes_factors_smart{5}(:, j));
    end
end

%% все расстояния между аффинитивностями
cost_functions = { (@(I, AC) norm(I - AC, 'fro'));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 1));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0.5));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 2));
                   (@(I, AC) nmf_alpha_divergence(I, AC, -1)); };
%I_hat = A_smart_full * C_smart_full;
probesets = unique(inten(:, end));
costs_all = zeros(length(probesets), 6);
for i = 1:6
    for j = 1:length(probesets)
        fprintf('%d %d\n', i, j);
        costs_all(j, i) = cost_functions{i}(arrays_factors_smart{3}(inten(:, end) == probesets(j)), ...
            arrays_factors_smart{6}(inten(:, end) == probesets(j)));
    end
end

%% одна оптимизация, все пары критериев
cost_functions = { (@(I, AC) norm(I - AC, 'fro'));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 1));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 0.5));
                   (@(I, AC) nmf_alpha_divergence(I, AC, 2));
                   (@(I, AC) nmf_alpha_divergence(I, AC, -1)); };
costs_all = zeros(100, 6);
A = arrays_factors_smart{1};
I = inten(:, setdiff(1:100, partition_arrays));
C = max(eps, (A' * A) \ (A' * I));

%LearningRate_C = repmat(0.1 ./ sum(A, 1)', [1 size(C, 2)]);
%T = A' * log(bsxfun(@rdivide, I, A * C + 1e-10));
%T = exp(bsxfun(@times, LearningRate_C, T));
%C = bsxfun(@times, C, T);

I_hat = A * C;

for i = 1:6
    for j = 1:50
        costs_all(j, i) = cost_functions{i}(I(:, j), I_hat(:, j));
    end
end

%% все оптимизации, один критерий
C_1 = { probes_factors_fro{2}, ...
        probes_factors_smart{2}, ...
        probes_factors_KL{2}, ...
        probes_factors_hellinger{2}, ...
        probes_factors_chi_squared{2}};
C_2 = { probes_factors_fro{5}, ...
        probes_factors_smart{5}, ...
        probes_factors_KL{5}, ...
        probes_factors_hellinger{5}, ...
        probes_factors_chi_squared{5}};

costs_all = zeros(100, 5);
for i = 1:5
    for j = 1:100
        costs_all(j, i) = nmf_alpha_divergence(C_1{i}(:, j), C_2{i}(:, j), 0);
    end
end

%% зависимость ошибки от размера обучения
for sample_size = 10:10:100
    subsample = randsample(100, sample_size, false);
    inten = inten_full(:, [subsample 201:202]);
    [A, C, Avect] = calibrate_model(inten, factorize_smart);
    C_test = max(eps, (A' * A) \ (A' * inten_test(:, 1:100)));
	train_error = nmf_alpha_divergence(inten(:, 1:sample_size), A * C, 0);
    test_error = nmf_alpha_divergence(inten_test(:, 1:100), A * C_test, 0);
    save(['sample_size_' int2str(sample_size)], 'A', 'C', 'Avect', 'C_test', 'train_error', 'test_error');
end

%%
plot(10:10:100, train_error_all ./ (10:10:100)', 'LineWidth', 2);
hold on
plot(10:10:100, test_error_all / 100, 'r', 'LineWidth', 2);
legend('Train error', 'Test error', 'Location', 'Best')