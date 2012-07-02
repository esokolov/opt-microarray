function C = nmf_alpha_beta_fixedA(I, A, alpha, beta, maxIterCnt, eps)
    C = nmf_init_als_fixedA(I, A, eps);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            Z = bsxfun(@times, bsxfun(@power, I + eps, alpha), ...
                               bsxfun(@power, A * C + eps, beta - 1));
            T = bsxfun(@rdivide, A' * Z, A' * bsxfun(@power, A * C + eps, alpha + beta - 1) + eps);
            C = bsxfun(@times, C, bsxfun(@power, T + eps, 1 / alpha));
        else
            %error('No alpha = 0 please');
            T = A' * bsxfun(@times, bsxfun(@power, A * C + eps, alpha + beta - 1), ...
                log(bsxfun(@rdivide, I + eps, A * C + eps) + eps));
            T = bsxfun(@rdivide, T + eps, A' * bsxfun(@power, A * C + eps, alpha + beta - 1));
            C = bsxfun(@times, C, exp(T));
        end
        
        currQuality = nmf_alpha_beta_divergence(I, A * C, alpha, beta);
        if (nmf_check_stopping_criteria(I, A, C, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end