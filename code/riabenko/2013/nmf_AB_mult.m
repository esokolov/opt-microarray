function [A X isConverged] = nmf_AB_mult(P, k, alpha, beta, maxIterCnt, eps)    
    eps_nnz = 1e-12;
    minIterCnt = 50;    

    [A X] = nmf_init_svd(P, k, eps_nnz);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha ~= 0)
            Z = bsxfun(@times, bsxfun(@power, P + eps_nnz, alpha), ...
                               bsxfun(@power, A * X + eps_nnz, beta - 1));
            T = bsxfun(@rdivide, A' * Z, A' * bsxfun(@power, A * X + eps_nnz, alpha + beta - 1) + eps_nnz);
            X = bsxfun(@times, X, bsxfun(@power, T + eps_nnz, 1 / alpha));
            
            Z = bsxfun(@times, bsxfun(@power, P + eps_nnz, alpha), ...
                               bsxfun(@power, A * X + eps_nnz, beta - 1));
            T = bsxfun(@rdivide, Z * X', bsxfun(@power, A * X + eps_nnz, alpha + beta - 1) * X' + eps_nnz);
            A = bsxfun(@times, A, bsxfun(@power, T + eps_nnz, 1 / alpha));
        else
            %error('No alpha = 0 please');
            T = A' * bsxfun(@times, bsxfun(@power, A * X + eps_nnz, alpha + beta - 1), ...
                log(bsxfun(@rdivide, P, A * X + eps_nnz) + eps_nnz));
            T = bsxfun(@rdivide, T + eps_nnz, A' * bsxfun(@power, A * X + eps_nnz, alpha + beta - 1));
            X = bsxfun(@times, X, exp(T));
            
            T = bsxfun(@times, bsxfun(@power, A * X + eps_nnz, alpha + beta - 1), ...
                log(bsxfun(@rdivide, P, A * X + eps_nnz) + eps_nnz)) * X';
            T = bsxfun(@rdivide, T + eps_nnz, bsxfun(@power, A * X + eps_nnz, alpha + beta - 1) * X');
            A = bsxfun(@times, A, exp(T));
        end
        
        [A X] = nmf_normalize_l2(A, X);
        
        currQuality = ABdiv(P, A * X, alpha, beta);
        if (currIter>minIterCnt && nmf_stopping(P, A, X, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        %вывести ещё score
        fprintf('%d: %f\n', currIter, currQuality);
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    X(isnan(X)) = 0;
end