function [A C] = nmf_smart(I, G, maxIterCnt, eps, LearningRate_A, LearningRate_C)
    [A C] = nmf_init_als(I, G, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        LearningRate_A = repmat(1 ./ sum(C, 2)', [size(A, 1) 1]);
        LearningRate_C = repmat(1 ./ sum(A, 1)', [1 size(C, 2)]);
        
        T = A' * log(bsxfun(@rdivide, I, A * C + eps));
        T = exp(bsxfun(@times, LearningRate_C, T));
        C = bsxfun(@times, C, T);
        
        T = log(bsxfun(@rdivide, I, A * C + eps)) * C';
        T = exp(bsxfun(@times, LearningRate_A, T));
        A = bsxfun(@times, A, T);
        
        [A C] = nmf_normalize(A, C);
        
        currQuality = nmf_alpha_divergence(I, A * C, 0);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end