function [A C isConverged] = nmf_alpha_beta_weighted(I, W, G, alpha, beta, maxIterCnt, eps)
    eps_nnz = 1e-12;
    minIterCnt = 50;    

    [A C] = nmf_init_als(I, G, eps_nnz);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            num = A' * bsxfun(@times, bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), bsxfun(@power, A * C + eps_nnz, beta - 1)), W);
            denom = A' * bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W);
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            C = bsxfun(@times, C, bsxfun(@power, T + eps_nnz, 1 / alpha));
            
            num = bsxfun(@times, bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), bsxfun(@power, A * C + eps_nnz, beta - 1)), W) * C';
            denom = bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W) * C';
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            A = bsxfun(@times, A, bsxfun(@power, T + eps_nnz, 1 / alpha));
        else
            %error('No alpha = 0 please');
            num = A' * bsxfun(@times, bsxfun(@times, log(bsxfun(@rdivide, I, A * C + eps_nnz) + eps_nnz), bsxfun(@power, A * C + eps_nnz, beta - 1)), W);
            denom = A' * bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W);
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            C = bsxfun(@times, C, exp(T));
            
            num = bsxfun(@times, bsxfun(@times, log(bsxfun(@rdivide, I, A * C + eps_nnz) + eps_nnz), bsxfun(@power, A * C + eps_nnz, beta - 1)), W) * C';
            denom = bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), W) * C';
            T = bsxfun(@rdivide, num, denom + eps_nnz);
            A = bsxfun(@times, A, exp(T));
        end
        
        [A C] = nmf_normalize_prod(A, C);
        
        currQuality = nmf_alpha_beta_divergence(I.*W, (A * C).*W, alpha, beta);
        if (currIter>minIterCnt && nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        %fprintf('%d: %f %e\n', currIter, currQuality, C(1689));
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end