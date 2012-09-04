function [dist_A, dist_B] = ...
        nonlinear_get_AB_error(inten, inten_sliced, inten_genes_idx, ...
        partition_arrays, ...
        alpha, beta, maxIterCnt, eps, reg_best)
    
    sampleSize = size(inten, 2) - 2;
    
    % разбиение по чипам
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

    tmp = abs(Avect1_arrays - Avect2_arrays) ./ (Avect1_arrays + Avect2_arrays);
    %tmp = tmp(~isnan(tmp));
    %tmp = tmp(~isinf(tmp));
    dist_A = sum(tmp>0.5)+sum(isnan(tmp))+sum(isinf(tmp));
    
    tmp = abs(Bvect1_arrays - Bvect2_arrays) ./ (Bvect1_arrays + Bvect2_arrays);
    %tmp = tmp(~isnan(tmp));
    %tmp = tmp(~isinf(tmp));
    dist_B = sum(tmp>0.5)+sum(isnan(tmp))+sum(isinf(tmp));
end