function [qual_A, outliers_A, qual_B, outliers_B] = ...
        get_quality_cv_frankenstein(inten, inten_sliced, inten_genes_idx, inten_test_sliced, ...
        partition_arrays, partition_probes, partition_probes_compl, partition_probes_sliced, partition_probes_sliced_compl, ...
        alpha, beta, maxIterCnt, eps, reg_best)
    
    sampleSize = size(inten, 2) - 2;
    
    % ��������� �� �����
    sample1 = inten(:, [partition_arrays (size(inten, 2)-1:size(inten, 2))]);
    sample2 = inten(:, [setdiff(1:sampleSize, partition_arrays) (size(inten, 2)-1:end)]);
    sample1_sliced = inten_sliced;
    sample2_sliced = inten_sliced;
    for i = 1:length(inten_sliced)
        sample1_sliced{i} = sample1_sliced{i}(:, partition_arrays);
        sample2_sliced{i} = sample2_sliced{i}(:, setdiff(1:sampleSize, partition_arrays));
    end

    fprintf(' 1');
    [~, ~, ~, Avect1_arrays, Bvect1_arrays] = nonlinear_calibrate_model_short(sample1, sample1_sliced, inten_genes_idx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, reg_best, 1));
    fprintf(' 2');
    [~, ~, ~, Avect2_arrays, Bvect2_arrays] = nonlinear_calibrate_model_short(sample2, sample2_sliced, inten_genes_idx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, reg_best, 1));

%     tmp = Avect1_arrays - Avect2_arrays;
%     tmp = tmp(~isnan(tmp));
%     tmp = tmp(~isinf(tmp));
%     mad_A = mad_my(tmp, 1);
%     outliers_A = sum(abs(tmp) > 3 * std(tmp));
    tmp = abs(Avect1_arrays - Avect2_arrays) ./ (Avect1_arrays + Avect2_arrays);
    tmp = tmp(~isnan(tmp));
    tmp = tmp(~isinf(tmp));
    qual_A = sum(tmp > 0.5);
    outliers_A = sum(abs(tmp) > 3 * std(tmp));
    
%     tmp = Bvect1_arrays - Bvect2_arrays;
%     tmp = tmp(~isnan(tmp));
%     tmp = tmp(~isinf(tmp));
%     mad_B = mad_my(tmp, 1);
%     outliers_B = sum(abs(tmp) > 3 * std(tmp));    
    tmp = abs(Bvect1_arrays - Bvect2_arrays) ./ (Bvect1_arrays + Bvect2_arrays);
    tmp = tmp(~isnan(tmp));
    tmp = tmp(~isinf(tmp));
    qual_B = sum(tmp > 0.5);
    outliers_B = sum(abs(tmp) > 3 * std(tmp));
    
%     % ��������� �� ������
%     sample1 = inten(partition_probes, :);
%     sample2 = inten(partition_probes_compl, :);
%     sample1_sliced = inten_sliced;
%     sample2_sliced = inten_sliced;
%     sample1_genes_idx = partition_probes_sliced;
%     sample2_genes_idx = partition_probes_sliced;
%     sample1_last_idx = 0;
%     sample2_last_idx = 0;
%     for i = 1:length(inten_sliced)
%         if (isempty(partition_probes_sliced{i}))
%             continue;
%         end
%         sample1_sliced{i} = sample1_sliced{i}(partition_probes_sliced{i}, :);
%         sample2_sliced{i} = sample2_sliced{i}(partition_probes_sliced_compl{i}, :);
%         sample1_genes_idx{i} = (sample1_last_idx + 1):(sample1_last_idx + length(partition_probes_sliced{i}));
%         sample2_genes_idx{i} = (sample2_last_idx + 1):(sample2_last_idx + length(partition_probes_sliced_compl{i}));
%         sample1_last_idx = sample1_last_idx + length(partition_probes_sliced{i});
%         sample2_last_idx = sample2_last_idx + length(partition_probes_sliced_compl{i});
%     end    
%     
%     for i = length(inten_sliced):-1:1
%         if (isempty(partition_probes_sliced{i}))
%             sample1_sliced(i) = [];
%             sample2_sliced(i) = [];
%             sample1_genes_idx(i) = [];
%             sample2_genes_idx(i) = [];
%         end
%     end
%     
%     fprintf(' 3');
%     [A1_probes, B1_probes, C1_probes, ~, A1_slices_probes, B1_sliced_probes, ~, ~, Avect1_probes, Bvect1_probes] = ...
%         nonlinear_calibrate_frankenstein_model(sample1, sample1_sliced, inten_test_sliced, inten_genes_idx, ...
%         alpha, beta, maxIterCnt, eps);
%     fprintf(' 4');
%     [A2_probes, B2_probes, C2_probes, ~, A2_slices_probes, B2_sliced_probes, ~, ~, Avect2_probes, Bvect2_probes] = ...
%         nonlinear_calibrate_frankenstein_model(sample2, sample2_sliced, inten_test_sliced, inten_genes_idx, ...
%         alpha, beta, maxIterCnt, eps);
% 
% 
%     tmp = C1_probes(:) - C2_probes(:);
%     tmp = tmp(~isnan(tmp));
%     tmp = tmp(~isinf(tmp));
%     mad_C = mad_my(tmp, 1);
%     outliers_C = sum(abs(tmp) > 3 * std(tmp));
%     
    
end