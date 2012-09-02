function [A B C isConverged] = nonlinear_alpha_beta(I, alpha, beta, maxIterCnt, eps, alpha_C, alpha_B, use_term_criteria)
    if (nargin < 8)
        use_term_criteria = true;
    end
        
    %[A B C] = nonlinear_init_als(I, eps);
    
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    
    minIterCnt = 50;
    eps_nnz = 1e-12;
    %[A C] = nmf_normalize_prod(A, C);
    
    %A = rand(size(I, 1), 1);
    %C = rand(1, size(I, 2));
    %B = rand(size(I, 1), 1);
    
    B = zeros(size(A)) + eps_nnz;
    
    isConverged = 1;
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha == 0 && beta == 0)
            %F = (A * C) ./ (1 + B * C);
            A = exp((1 / size(I, 2)) * sum(log(I + eps_nnz) - log(bsxfun(@rdivide, C, 1 + B * C) + eps_nnz), 2));
            A = max(A, 0);
            
            F = (A * C) ./ (1 + B * C);
            C = C + 0.5 * (sum((1 ./ (eps_nnz + bsxfun(@plus, C, B * (C .^ 2)))) .* (log(I + eps_nnz) - log(F + eps_nnz)), 1) + alpha_C * C) ./ ...
                (sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz)) + 1), 1) + alpha_C/2);            
            C = max(C, 0);
            
            F = (A * C) ./ (1 + B * C);
            B = B - 0.5 * (sum(bsxfun(@rdivide, C, 1 + B * C) .* (log(I + eps_nnz) - log(F + eps_nnz)), 2) + alpha_B * B) ./ ...
                (sum(bsxfun(@times, C .^ 2, (1 + B * C) .^ 2) .* (log(F + eps_nnz) - log(I + eps_nnz) + 1), 2) + alpha_B);
            B = max(B, 0);
        elseif (alpha == 0)
            F = (A * C) ./ (1 + B * C);
            A = A - 0.5 * ((1 ./ (A + eps_nnz)) .* sum(((F + eps_nnz) .^ beta) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                ((1 ./ (eps_nnz + A .^ 2)) .* sum(((F + eps_nnz) .^ beta) .* ((beta - 1) * log((F + eps_nnz) ./ (I + eps_nnz)) + 1), 2));
            A = max(A, 0);
            
            F = (A * C) ./ (1 + B * C);
            C = C - 0.5 * (sum(((F + eps_nnz) .^ (beta + 1)) .* log((F + eps_nnz) ./ (I + eps_nnz)) .* (1 ./ (eps_nnz + A * (C .^ 2))), 1) + alpha_C * C) ./ ...
                (sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((F + eps_nnz) .^ beta) .* ...
                (1 - (2 * B * C - beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz))), 1) + alpha_C/2);
            C = max(C, 0);
            
            F = (A * C) ./ (1 + B * C);
            B = B + 0.5 * (alpha_B * B + (1 ./ (A + eps_nnz)) .* sum(((F + eps_nnz) .^ (beta + 1)) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (sum(bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2) .* ((F + eps_nnz) .^ beta) .* ...
                ((beta + 1) * log((F + eps_nnz) ./ (I + eps_nnz)) + 1), 2) + alpha_B/2);
            B = max(B, 0);
        elseif (beta == 0)
            A = (sum((I + eps_nnz) .^ alpha, 2) ./ sum((eps_nnz + bsxfun(@rdivide, C, 1 + B * C)) .^ alpha, 2)) .^ (1 / alpha);
            A = max(A, 0);
            
            C = C - 0.5 * (sum((((eps_nnz + (A * C) ./ (1 + B * C)) .^ alpha) - (eps_nnz + I) .^ alpha) ./ ...
                bsxfun(@times, C, 1 + B * C), 1) + alpha_C * C) ./ ...
                (sum(((alpha - 1 - 2 * B * C) .* ((eps_nnz + (A * C) ./ (1 + B * C)) .^ alpha) + ...
                (1 + 2 * B * C) .* ((I + eps_nnz) .^ alpha)) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + alpha_C/2);
            C = max(C, 0);
            
            B = B - 0.5 * (sum(bsxfun(@times, C, (I + eps_nnz) .^ alpha - (eps_nnz + (A * C) ./ (1 + B * C)) .^ alpha) ./ ...
                (1 + B * C), 2) + alpha_B * B) ./ ...
                (sum(bsxfun(@times, C .^ 2, (alpha + 1) * (eps_nnz + (A * C) ./ (1 + B * C)) .^ alpha - (I + eps_nnz) .^ alpha) ./ ...
                ((1 + B * C) .^ 2), 2) + alpha_B/2);
            B = max(B, 0);
        elseif (alpha == -beta)
            A = ((1 / size(I, 2)) * sum((bsxfun(@rdivide, I .* (1 + B * C), C) + eps_nnz) .^ alpha, 2)) .^ (1 / alpha);
            A = max(A, 0);
            
            C = C - 0.5 * (sum((1 - (eps_nnz + (I .* (1 + B * C)) ./ (A * C)) .^ alpha) ./ (eps_nnz + bsxfun(@times, C, 1 + B * C)), 1) + alpha_C * C ./ quantile(I, 0.9))./ ...
                (sum(((((I .* (1 + B * C)) ./ (A * C) + eps_nnz) .^ alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + alpha_C ./ (2 * quantile(I, 0.9)));
            C = max(C, 0);
            
            B = B - 0.5 * (sum(((eps_nnz + bsxfun(@rdivide, I, A)) .^ alpha) .* ((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)) .^ (alpha - 1)) - ...
                bsxfun(@rdivide, C, 1 + B * C), 2) + alpha_B * B ./ quantile(I', 0.9)') ./ ...
                (sum((alpha - 1) * (((eps_nnz + bsxfun(@rdivide, I, A)) .^ alpha) .* ((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)) .^ (alpha - 2))) + ...
                bsxfun(@rdivide, C, 1 + B * C) .^ 2, 2) + alpha_B ./ (2 * quantile(I', 0.9)'));
            B = max(B, 0);
        else
            % optimizing A
            F = bsxfun(@rdivide, C, 1 + B * C);
            A = (sum(((I + eps_nnz) .^ alpha) .* ((F + eps_nnz) .^ beta), 2) ./ sum((F + eps_nnz) .^ (alpha + beta), 2)) .^ (1 / alpha);
            A = max(A, 0);

            % optimizing C
            F = (A * C) ./ (1 + B * C);
            C = C - 0.5 * (sum((1 ./ (eps_nnz + A * (C .^ 2))) .* ((F + eps_nnz) .^ (beta + 1)) .* ((I + eps_nnz) .^ alpha - (F + eps_nnz) .^ alpha), 1) + alpha_C * C) ./ ...
                (sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((F + eps_nnz) .^ beta) .* (((F + eps_nnz) .^ alpha) .* (2 * B * C - alpha - beta + 1) - ...
                ((I + eps_nnz) .^ alpha) .* (2 * B * C - beta + 1)), 1) + alpha_C/2);
            C = max(C, 0);

            % optimizing B
            F = (A * C) ./ (1 + B * C);
            B = B + 0.5 * ((1 ./ (A + eps_nnz)) .* sum(((F + eps_nnz) .^ (beta + 1)) .* ((I + eps_nnz) .^ alpha - (F + eps_nnz) .^ alpha), 2) + alpha_B * B) ./ ...
                (sum(bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2) .* ((F + eps_nnz) .^ beta) .* ((beta + 1) * ((I + eps_nnz) .^ alpha) - ...
                (alpha + beta + 1) * ((F + eps_nnz) .^ alpha)), 2) + alpha_B/2);
            B = max(B, 0);
        end
             
        
        % routines
        if (alpha_B == 0 && alpha_C == 0)
            [A, B, C] = nonlinear_normalize_prod(A, B, C);
        end
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        
        %if (currIter > 1 && currQuality > prevQuality)
        %    isConverged = false;
        %    break;
        %end
        
        %if (currQuality > prevQuality && sum(C == 0) > 0)
        %    fprintf('%d\n', find(C == 0));
        %    I = I(:, setdiff(1:size(I, 2), find(C == 0)));
        %    [A B C] = nonlinear_alpha_beta(I, alpha, beta, maxIterCnt, eps, opt_method_C, opt_method_B, use_term_criteria);
        %    break;
        %end
        
        if (currIter>minIterCnt && use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        %if (currIter > 10000)
        %    break;
        %end
        fprintf('%d: %f\n', currIter, currQuality);
        %fprintf('%d: %e\n', currIter, C(525));
        
    end
    
    %isConverged = (currIter < maxIterCnt);
    
    [A, B, C] = nonlinear_normalize_prod(A, B, C);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end
