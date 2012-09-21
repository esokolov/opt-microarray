function [C isConverged] = nmf_alpha_beta_fixedA_weighted(I, W, A, alpha, beta, maxIterCnt, eps)
    C = nmf_init_als_fixedA(I, A, eps);
    
    eps_nnz = 1e-12;
    minIterCnt = 50;
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            num = A' * bsxfun(@times, bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), bsxfun(@power, A * C + eps_nnz, beta - 1)), W);
            denom = A' * bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W);
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            C = bsxfun(@times, C, bsxfun(@power, T + eps_nnz, 1 / alpha));
        else
            %error('No alpha = 0 please');
            num = A' * bsxfun(@times, bsxfun(@times, log(bsxfun(@rdivide, I, A * C + eps_nnz) + eps_nnz), bsxfun(@power, A * C + eps_nnz, beta - 1)), W);
            denom = A' * bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W);
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            C = bsxfun(@times, C, exp(T));
        end
        
        currQuality = nmf_alpha_beta_divergence(I, A * C, alpha, beta);
        if (currIter>minIterCnt && nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end