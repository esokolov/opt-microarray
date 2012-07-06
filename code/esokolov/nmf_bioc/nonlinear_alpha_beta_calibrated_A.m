%function [A B C isConverged] = nonlinear_alpha_beta(I, alpha, beta, maxIterCnt, eps)
function [C isConverged] = nonlinear_alpha_beta_calibrated_A(I, A_calibrated, alpha, beta, maxIterCnt, eps)
    %[A B C] = nonlinear_init_als(I, eps);
    %[A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    %[A C] = nmf_normalize_prod(A, C);
    C = nmf_alpha_beta_fixedA(I, A_calibrated, alpha, beta, maxIterCnt, eps);
    
    A = A_calibrated;
    B = zeros(size(A));
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        % optimizing A
        A = size(I, 2) ./ sum(bsxfun(@rdivide, C, I .* (1 + B * C)), 2);
        %F = (A * C) ./ (1 + B * C);
        %D = bsxfun(@rdivide, C, 1 + B * C);
        %A = A .* (sum((D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1))), 2) ./ ...
        %    sum(eps + ((F + eps) .^ (alpha + beta - 1)) .* D, 2)) .^ (1 / alpha);
         
        % optimizing C
        %C = projected_grad(C, maxIterCnt, 1e-15, ...
        %    @(C_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B, C_arg), alpha, beta), ...
        %    @(C_arg) nonlinear_alpha_beta_grad_C_special(I, A, B, C_arg));
        C = sum(1 ./ (1 + B * C), 1) ./ sum(bsxfun(@rdivide, A, I .* ((1 + B * C) .^ 2)), 1);
        
        %F = (A * C) ./ (1 + B * C);
        %D = bsxfun(@rdivide, A, (1 + B * C) .^ 2);
        %C = C .* (sum(D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)), 1) ./ ...
        %    sum(eps + D .* ((F + eps) .^ (alpha + beta - 1)), 1)) .^ (1 / alpha);
        
        % optimizing B
        B = projected_grad(B, maxIterCnt, 1e-15, ...
            @(B_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B_arg, C), alpha, beta), ...
            @(B_arg) nonlinear_alpha_beta_grad_B_special(I, A, B_arg, C));
        
%         B = projected_grad(B, maxIterCnt, 1e-15, ...
%             @(B_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B_arg, C), alpha, beta), ...
%             @(B_arg) nonlinear_alpha_beta_grad_B(I, A, B_arg, C, alpha, beta));
        %B = projected_grad_B(I, A, B, C, alpha, beta, maxIterCnt, 1e-15);        
        
        % routines
        [A, B, C] = nonlinear_normalize_prod(A, B, C);
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        if (nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        
        %fprintf('%d: %f\n', currIter, currQuality);
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function X_new = projected_grad(X, maxIterCnt, eps, cost_func, grad_func)
    sigma = 1/4;
    eta_dec = 0.5;
    eta_init = max(1000 * max(X), 1);
    
    prevCost = -1;
    for currIter = 1:maxIterCnt
        currCost = cost_func(X);
        currGrad = grad_func(X);
        
        if (sum(sum(abs(currGrad))) < eps)
            break;
        end
        if (currIter > 1 && prevCost - currCost < eps)
            break;
        end
        prevCost = currCost;
        
        eta = eta_init;
        while (eta > eps)
            X_new = max(X - eta * currGrad, 0);
            
            newCost = cost_func(X_new);
            
            armijo_cond = ((newCost - currCost) <= sigma * (currGrad(:)' * (X_new(:) - X(:))));
            if armijo_cond
                break;
            else
                eta = eta * eta_dec;
            end
        end
        
        X = X_new;
        
        break;
    end
    
    X_new = X;
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
        break;
    end
    
    B_new = B;
end

function grad = nonlinear_alpha_beta_grad_B_special(I, A, B, C)
    grad = -sum( ((A * C) ./ I - (1 + B * C)) .* bsxfun(@rdivide, C, (1 + B * C) .^ 2), 2);
end

function grad = nonlinear_alpha_beta_grad_C_special(I, A, B, C)
    grad = sum( (bsxfun(@rdivide, A, I) - bsxfun(@rdivide, 1 + B * C, C)) .* bsxfun(@rdivide, 1, (1 + B * C) .^ 2), 1);
end

% function grad = nonlinear_alpha_beta_grad_B(I, A, B, C, alpha, beta)
%     F = (A * C) ./ (1 + B * C);
%     D = (A * (C .^ 2)) ./ ((1 + B * C) .^ 2);
%     grad = sum(-(((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)) .* D) + ...
%         ((F + eps) .^ (alpha + beta - 1) .* D), 2);
% end