% opt_method:
%   mult
%   quadratic
%   projected_grad
%   equation_solver
%   sqp
%   trust-region-reflective
%   active-set
%   interior-point
function [A B C isConverged Quality] = nonlinear_alpha_beta_wq(I, alpha, beta, maxIterCnt, eps, opt_method_C, opt_method_B, use_term_criteria)
    if (nargin < 6)
        opt_method_C = 'mult';
        opt_method_B = 'mult';
    elseif (nargin < 7)
        opt_method_B = 'mult';
    elseif (nargin < 8)
        use_term_criteria = true;
    end
        
    %[A B C] = nonlinear_init_als(I, eps);
    
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    
    %[A C] = nmf_normalize_prod(A, C);
    
    %A = rand(size(I, 1), 1);
    %C = rand(1, size(I, 2));
    %B = rand(size(I, 1), 1);
    
    B = zeros(size(A)) + eps;
    
    Quality = zeros(maxIterCnt,1);
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        % optimizing A
        A = size(I, 2) ./ sum(bsxfun(@rdivide, C, I .* (1 + B * C)), 2);
        %F = (A * C) ./ (1 + B * C);
        %D = bsxfun(@rdivide, C, 1 + B * C);
        %A = A .* (sum((D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1))), 2) ./ ...
        %    sum(eps + ((F + eps) .^ (alpha + beta - 1)) .* D, 2)) .^ (1 / alpha);
         
        % optimizing C
        if strcmp(opt_method_C, 'projected_grad')
            C = projected_grad(C, maxIterCnt, 1e-15, ...
                @(C_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B, C_arg), alpha, beta), ...
                @(C_arg) nonlinear_alpha_beta_grad_C_special(I, A, B, C_arg));
        elseif strcmp(opt_method_C, 'mult')
            C = sum(1 ./ (1 + B * C), 1) ./ sum(bsxfun(@rdivide, A, I .* ((1 + B * C) .^ 2)), 1);
        elseif strcmp(opt_method_C, 'quadratic')
            C = C - 0.5 * sum(bsxfun(@rdivide, A, I .* ((1 + B * C) .^ 2)) - 1 ./ bsxfun(@plus, C, B * (C .^ 2)), 1) ./ ...
                sum((B * C + 0.5) ./ (bsxfun(@times, C, 1 + B * C) .^ 2) - bsxfun(@rdivide, A .* B, I .* ((1 + B * C) .^ 3)), 1);
            C = max(C, 0);
        elseif strcmp(opt_method_C, 'equation_solver')
            for i = 1:length(C)
                f = @(c) (-sum(1 ./ (c .* (1 + B * c))) + sum(A ./ (I(:, i) .* (1 + B * C(i)))));
                C(i) = max(0, fsolve(f, C(i), optimset('Display', 'off')));
            end
        else        
            for i = 1:length(C)
                f_cost = @(c) (sum((A .* C(i)) ./ (I(:, i) .* ((1 + B * C(i)) .^ 2))) - sum(log((A .* c) ./ (1 + B * c))));
                f_grad = @(c) (-sum(1 ./ (c .* (1 + B * c))) + sum(A ./ (I(:, i) .* (1 + B * C(i)))));
                f = @(c) fff(c, f_cost, f_grad);
                C(i) = fmincon(f, C(i), [], [], [], [], 0, [], [], ...
                    optimset('Algorithm', opt_method_C, 'GradObj', 'on', 'Display', 'off'));
                %C(i) = C_new + 0.5 * (C_new - C(i));
            end
        end
        
%         F = (A * C) ./ (1 + B * C);
%         D = bsxfun(@rdivide, A, (1 + B * C) .^ 2);
%         C = C .* (sum(D .* ((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)), 1) ./ ...
%             sum(eps + D .* ((F + eps) .^ (alpha + beta - 1)), 1)) .^ (1 / alpha);
        
        % optimizing B
        if strcmp(opt_method_B, 'projected_grad')
            B = projected_grad(B, maxIterCnt, 1e-15, ...
                @(B_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B_arg, C), alpha, beta), ...
                @(B_arg) nonlinear_alpha_beta_grad_B_special(I, A, B_arg, C));
        elseif strcmp(opt_method_B, 'mult')
            B = B .* (sum((A * (C .^ 2)) ./ (I .* ((1 + B * C) .^ 2)), 2) ./ ...
                sum(bsxfun(@rdivide, C, 1 + B * C), 2));
        elseif strcmp(opt_method_B, 'quadratic')
            rpk = A * C ./ I;
            B = B - 0.5 * sum(bsxfun(@rdivide, C, 1 + B * C) - bsxfun(@times, rpk, C) ./ ((1 + B * C) .^ 2), 2) ./ ...
                sum(bsxfun(@times, rpk, C .^ 2) ./ ((1 + B * C) .^ 3) - 0.5 * bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2), 2);
            B = max(B, 0);
        elseif strcmp(opt_method_C, 'equation_solver')
            for i = 1:length(B)
                f = @(b) (sum((A(i) * (C .^ 2)) ./ (I(i, :) .* ((1 + b * C) .^ 2))) - sum(C ./ (1 + B(i) * C)));
                B(i) = max(0, fzero(f, B(i), optimset('Display', 'off')));
            end
        else
            for i = 1:length(B)
                f_cost = @(b) (sum((C * b) ./ (1 + B(i) * C)) + sum((A(i) * C) ./ (I(i, :) .* (1 + b * C))));
                f_grad = @(b) (-sum((A(i) * (C .^ 2)) ./ (I(i, :) .* ((1 + b * C) .^ 2))) + sum(C ./ (1 + B(i) * C)));
                f = @(b) fff(b, f_cost, f_grad);
                B(i) = fmincon(f, B(i), [], [], [], [], 0, [], [], ...
                    optimset('Algorithm', 'sqp', 'GradObj', 'on', 'Display', 'off'));
            end
        end
        
        %B = B .* (sum((A * (C .^ 2)) ./ (I .* ((1 + B * C) .^ 2)), 2) ./ ...
        %    sum(bsxfun(@rdivide, C, 1 + B * C), 2));
        
        % quadratic approximation
        
        
        
        
%         B = projected_grad(B, maxIterCnt, 1e-15, ...
%             @(B_arg) nmf_alpha_beta_divergence(I, langmuir_func(A, B_arg, C), alpha, beta), ...
%             @(B_arg) nonlinear_alpha_beta_grad_B(I, A, B_arg, C, alpha, beta));
        %B = projected_grad_B(I, A, B, C, alpha, beta, maxIterCnt, 1e-15);        
        
        % routines
        [A, B, C] = nonlinear_normalize_prod(A, B, C);
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        if (use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        Quality(currIter) = currQuality;
        %if (currIter > 10000)
        %    break;
        %end
        fprintf('%d: %f\n', currIter, currQuality);
    end
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function [a1 a2] = fff(c, f_cost, f_grad)
    a1 = f_cost(c);
    a2 = f_grad(c);
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