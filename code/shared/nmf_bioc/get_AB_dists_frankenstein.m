dist_A = zeros(length(alpha_range), length(beta_range));
dist_B = zeros(length(alpha_range), length(beta_range));

for i = 1:length(alpha_range)
    for j = 1:length(beta_range)
        if (alpha_range(i) < -beta_range(j))
            continue;
        end
        
        fprintf('Alpha = %f, Beta = %f\n', alpha_range(i), beta_range(j));
        
        [dist_A(i,j), dist_B(i,j)] = nonlinear_get_AB_error(inten_train, inten_train_sliced, inten_full_idx, ...
        partition_arrays, alpha_range(i), beta_range(j), maxIterCnt, eps, reg_best(i, j));
    end
end