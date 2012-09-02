function [A C isConverged] = nmf_alpha_beta(I, G, alpha, beta, maxIterCnt, eps)    
    eps_nnz = 1e-12;
    minIterCnt = 50;    

    [A C] = nmf_init_als(I, G, eps_nnz);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            Z = bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), ...
                               bsxfun(@power, A * C + eps_nnz, beta - 1));
            T = bsxfun(@rdivide, A' * Z, A' * bsxfun(@power, A * C + eps_nnz, alpha + beta - 1) + eps_nnz);
            C = bsxfun(@times, C, bsxfun(@power, T + eps_nnz, 1 / alpha));
            
            Z = bsxfun(@times, bsxfun(@power, I + eps_nnz, alpha), ...
                               bsxfun(@power, A * C + eps_nnz, beta - 1));
            T = bsxfun(@rdivide, Z * C', bsxfun(@power, A * C + eps_nnz, alpha + beta - 1) * C' + eps_nnz);
            A = bsxfun(@times, A, bsxfun(@power, T + eps_nnz, 1 / alpha));
        else
            %error('No alpha = 0 please');
            T = A' * bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I, A * C + eps_nnz) + eps_nnz));
            T = bsxfun(@rdivide, T + eps_nnz, A' * bsxfun(@power, A * C + eps_nnz, alpha + beta - 1));
            C = bsxfun(@times, C, exp(T));
            
            T = bsxfun(@times, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I, A * C + eps_nnz) + eps_nnz)) * C';
            T = bsxfun(@rdivide, T + eps_nnz, bsxfun(@power, A * C + eps_nnz, alpha + beta - 1) * C');
            A = bsxfun(@times, A, exp(T));
        end
        
        [A C] = nmf_normalize_prod(A, C);
        
        currQuality = nmf_alpha_beta_divergence(I, A * C, alpha, beta);
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