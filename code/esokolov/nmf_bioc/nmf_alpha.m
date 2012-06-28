function [A C] = nmf_alpha(I, G, alpha, maxIterCnt, eps)
    [A C] = nmf_init_als(I, G, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        T = bsxfun(@rdivide, I + eps, A * C + eps);
        T = A' * bsxfun(@power, T, alpha);
        C = bsxfun(@times, C, bsxfun(@power, T, 1 / alpha));
        
        T = bsxfun(@rdivide, I + eps, A * C + eps);
        T = bsxfun(@power, T, alpha) * C';
        A = bsxfun(@times, A, bsxfun(@power, T, 1 / alpha));
        
        [A C] = nmf_normalize(A, C);
        
        currQuality = nmf_alpha_divergence(I, A * C, alpha);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end