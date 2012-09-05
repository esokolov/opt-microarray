function LOO_error = nmf_rma_get_LOO_error(inten_sliced, A_sliced, C_full)
    OUT_CNT = 5;

    LOO_error = 0;
    
    parfor i = 1:length(inten_sliced)
        LOO_error_curr = 0;
        
        if (size(inten_sliced{i}, 1) > OUT_CNT)
            out_idx = randsample(size(inten_sliced{i}, 1), OUT_CNT)';
        else
            out_idx = 1:size(inten_sliced{i}, 1);
        end
        
        C_diff = zeros(1, size(C_full, 2));
        
        for out_num = out_idx %1:size(inten_sliced{i}, 1)
            inten_curr = inten_sliced{i}(setdiff(1:size(inten_sliced{i}, 1), out_num), :);
            A_curr = A_sliced{i}(setdiff(1:size(inten_sliced{i}, 1), out_num));
            [A_curr] = nmf_normalize_prod(A_curr, 1);
            C_curr = nmf_rma_fixedA(inten_curr, A_curr);
            
            C_diff = C_diff + abs(C_curr - C_full(i, :)) ./ (C_curr + C_full(i, :));
            %LOO_error_curr = LOO_error_curr + median(abs(C_curr - C_full(i)));
        end
        
        C_diff = C_diff / length(out_idx);
        
        LOO_error_curr = mean(C_diff);
        
        LOO_error = LOO_error + LOO_error_curr;
    end
    
    LOO_error = LOO_error / length(inten_sliced);
end