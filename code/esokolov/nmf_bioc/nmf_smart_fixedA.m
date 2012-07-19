function C = nmf_smart_fixedA(I, A, maxIterCnt, eps)
    C = nmf_init_als_fixedA(I, A, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        LearningRate_C = repmat(0.1 ./ sum(A, 1)', [1 size(C, 2)]);
        
        T = A' * log(bsxfun(@rdivide, I + eps, A * C + eps));
        T = exp(bsxfun(@times, LearningRate_C, T));
        C = bsxfun(@times, C, T);
        
        currQuality = nmf_alpha_divergence(I, A * C, 0);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    C(isnan(C)) = 0;
end