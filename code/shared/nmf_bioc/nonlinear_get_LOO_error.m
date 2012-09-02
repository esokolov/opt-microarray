function LOO_error = nonlinear_get_LOO_error(inten_sliced, A_sliced, B_sliced, C_full, alpha_C, alpha, beta, maxIterCnt, eps)
    OUT_CNT = 5;

    LOO_error = 0;
    
    parfor i = 1:length(inten_sliced)
        LOO_error_curr = 0;
        
        if (size(inten_sliced{i}, 1) > OUT_CNT)
            out_idx = randsample(size(inten_sliced{i}, 1), OUT_CNT)';
        else
            out_idx = 1:size(inten_sliced{i}, 1);
        end
        
        for out_num = out_idx %1:size(inten_sliced{i}, 1)
            inten_curr = inten_sliced{i}(setdiff(1:size(inten_sliced{i}, 1), out_num), :);
            A_curr = A_sliced{i}(setdiff(1:size(inten_sliced{i}, 1), out_num));
            B_curr = B_sliced{i}(setdiff(1:size(inten_sliced{i}, 1), out_num));
            [A_curr, B_curr] = nonlinear_normalize_prod(A_curr, B_curr, 1);
            C_curr = nonlinear_alpha_beta_LO_reg_fixedAB(inten_curr, A_curr, B_curr, alpha, beta, maxIterCnt, eps, alpha_C, 1);
            
            LOO_error_curr = LOO_error_curr + median(abs(C_curr - C_full(i)));
        end
        
        LOO_error = LOO_error + LOO_error_curr / size(inten_sliced{i}, 1);
    end
    
    LOO_error = LOO_error / length(inten_sliced);
end