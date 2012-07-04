function [A C isConverged] = nmf_alpha_beta_full(I, G, alpha, beta, maxIterCnt, eps)
    [A C] = nmf_init_als(I, G, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            for i = 1:50
                Z = bsxfun(@times, bsxfun(@power, I + eps, alpha), ...
                                   bsxfun(@power, A * C + eps, beta - 1));
                T = bsxfun(@rdivide, A' * Z, A' * bsxfun(@power, A * C + eps, alpha + beta - 1) + eps);
                C = bsxfun(@times, C, bsxfun(@power, T + eps, 1 / alpha));
            end
            
            for i = 1:50
                Z = bsxfun(@times, bsxfun(@power, I + eps, alpha), ...
                                   bsxfun(@power, A * C + eps, beta - 1));
                T = bsxfun(@rdivide, Z * C', bsxfun(@power, A * C + eps, alpha + beta - 1) * C' + eps);
                A = bsxfun(@times, A, bsxfun(@power, T + eps, 1 / alpha));
            end
        else
            %error('No alpha = 0 please');
            T = A' * bsxfun(@times, bsxfun(@power, A * C + eps, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I + eps, A * C + eps) + eps));
            T = bsxfun(@rdivide, T + eps, A' * bsxfun(@power, A * C + eps, alpha + beta - 1));
            C = bsxfun(@times, C, exp(T));
            
            T = bsxfun(@times, bsxfun(@power, A * C + eps, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I + eps, A * C + eps) + eps)) * C';
            T = bsxfun(@rdivide, T + eps, bsxfun(@power, A * C + eps, alpha + beta - 1) * C');
            A = bsxfun(@times, A, exp(T));
        end
        
        [A C] = nmf_normalize_prod(A, C);
        
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