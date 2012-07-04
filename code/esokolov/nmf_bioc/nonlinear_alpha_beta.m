function [A B C isConverged] = nonlinear_alpha_beta(I, alpha, beta, maxIterCnt, eps)
    %[A B C] = nonlinear_init_als(I, eps);
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    %[A C] = nmf_normalize_prod(A, C);
    B = zeros(size(A));
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        % optimizing A
        %C_prime = C ./ (1 + B * C);
        %A = nmf_alpha_beta_fixedC(I, C_prime, alpha, beta, maxIterCnt, eps);
        F = (A * C) ./ (1 + B * C);
        D = bsxfun(@rdivide, C, 1 + B * C);
        A = A .* (sum((D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1))), 2) ./ ...
            sum(eps + ((F + eps) .^ (alpha + beta - 1)) .* D, 2)) .^ (1 / alpha);
         
        % optimizing C
        F = (A * C) ./ (1 + B * C);
        D = bsxfun(@rdivide, A, (1 + B * C) .^ 2);
        C = C .* (sum(D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)), 1) ./ ...
            sum(eps + D .* ((F + eps) .^ (alpha + beta - 1)), 1)) .^ (1 / alpha);
        
        % optimizing B
        B = projected_grad_B(I, A, B, C, alpha, beta, maxIterCnt, 1e-15);        
        
        % routines
        [A, B, C] = nonlinear_normalize_prod(A, B, C);
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        if (nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function B_new = projected_grad_B(I, A, B, C, alpha, beta, maxIterCnt, eps)
    sigma = 1/4;
    eta_dec = 0.5;
    eta_init = max(1000 * max(B), 1);
    
    prevCost = -1;
    for currIter = 1:maxIterCnt
        currCost = nmf_alpha_beta_divergence(I, langmuir_func(A, B, C), alpha, beta);
        currGrad = nonlinear_alpha_beta_grad_B(I, A, B, C, alpha, beta);
        
        if (sum(sum(abs(currGrad))) < eps)
            break;
        end
        if (currIter > 1 && prevCost - currCost < eps)
            break;
        end
        prevCost = currCost;
        
        eta = eta_init;
        while (eta > eps)
            B_new = max(B - eta * currGrad, 0);
            
            newCost = nmf_alpha_beta_divergence(I, langmuir_func(A, B_new, C), alpha, beta);
            
            armijo_cond = ((newCost - currCost) <= sigma * (currGrad(:)' * (B_new(:) - B(:))));
            if armijo_cond
                break;
            else
                eta = eta * eta_dec;
            end
        end
        
        B = B_new;
        %prevCost = newCost;
    end
    
    B_new = B;
end

function grad = nonlinear_alpha_beta_grad_B(I, A, B, C, alpha, beta)
    F = (A * C) ./ (1 + B * C);
    D = (A * (C .^ 2)) ./ ((1 + B * C) .^ 2);
    grad = sum(-(((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)) .* D) + ...
        ((F + eps) .^ (alpha + beta - 1) .* D), 2);
end