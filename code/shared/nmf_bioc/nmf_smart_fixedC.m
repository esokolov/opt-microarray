function A = nmf_smart_fixedC(I, C, maxIterCnt, eps)
    A = nmf_init_als_fixedC(I, C, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        LearningRate_A = repmat(1 ./ sum(C, 2)', [size(A, 1) 1]);
        
        T = log(bsxfun(@rdivide, I + eps, A * C + eps)) * C';
        T = exp(bsxfun(@times, LearningRate_A, T));
        A = bsxfun(@times, A, T);
        
        currQuality = nmf_alpha_divergence(I, A * C, 0);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    A(isnan(A)) = 0;
end