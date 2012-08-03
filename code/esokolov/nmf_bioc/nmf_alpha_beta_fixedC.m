function [A isConverged] = nmf_alpha_beta_fixedC(I, C, alpha, beta, maxIterCnt, eps)
    A = nmf_init_als_fixedC(I, C, eps);
    
    eps_nnz = 1e-12;
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            Z = bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), ...
                               bsxfun(@power, A * C + eps_nnz, beta - 1));
            T = bsxfun(@rdivide, Z * C', bsxfun(@power, A * C + eps_nnz, alpha + beta - 1) * C' + eps_nnz);
            A = bsxfun(@times, A, bsxfun(@power, T + eps_nnz, 1 / alpha));
        else
            T = bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I + eps_nnz, A * C + eps_nnz) + eps_nnz)) * C';
            T = bsxfun(@rdivide, T + eps_nnz, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1) * C');
            A = bsxfun(@times, A, exp(T));
        end
        
        currQuality = nmf_alpha_beta_divergence(I, A * C, alpha, beta);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end