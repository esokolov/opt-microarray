alpha_range = -2:0.25:2;
beta_range = -2:0.25:4;

mask = (norm_problems_cnt < 10) & (zero_influence_ratio < 1.1) & (outliers_influence_ratio < 1000) & (converged_cnt_full < 20);

quality_factors = cell(5, 1);
quality_factors{1} = concentration_dists .* (mask & (concentration_dists <= concentration_dists_mp));
quality_factors{2} = quality_arrays_mad .* (mask & (quality_arrays_mad <= quality_arrays_mad_mp));
quality_factors{3} = quality_probes_mad .* (mask & (quality_probes_mad <= quality_probes_mad_mp));
quality_factors{4} = overfitting_fro .* mask;
quality_factors{5} = overfitting_native .* mask;
quality_factors{6} = goodness_of_fit_fro .* mask;
quality_factors{7} = goodness_of_fit_geman .* mask;

weights = [1 1 0 1 1 0 0];

quality_final = zeros(size(quality_factors{1}));
for i = 1:length(quality_factors)
    quality_factors{i}(quality_factors{i} == 0) = Inf;
    quality_factors{i}(isnan(quality_factors{i})) = Inf;
    
    tmp = quality_factors{i}(:);
    [~, idx] = sort(tmp);
    tmp(idx) = 1:length(tmp);
    quality_factors{i} = reshape(tmp, size(quality_factors{i}));
    
    quality_final = quality_final + weights(i) * quality_factors{i};
end

[alpha_idx beta_idx] = ind2sub(size(quality_final), ...
    find(quality_final == min(min(quality_final)), 1, 'first'));

alpha = alpha_range(alpha_idx);
beta = alpha_range(beta_idx);

%%
for i = 1:7
    figure;
    imagesc(alpha_range, beta_range, quality_factors{i});
    hold on;
    set(gca, 'YDir', 'normal');
    plot([-2 2], [0 0]);
    plot([-3, 3], [4, -2]);
    switch i
        case 1
            title('MAD(D_c(C_1, C_2))');
        case 2
            title('MAD(A_1 - A_2)');
        case 3
            title('MAD(C_1 - C_2)');
        case 4
            title('Overfitting (frobenius norm)');
        case 5
            title('Overfitting (native divergence)');
        case 6
            title('Goodness of fit (frobenius norm)');
        case 7
            title('Goodness of fit (geman m-function)');
    end
end