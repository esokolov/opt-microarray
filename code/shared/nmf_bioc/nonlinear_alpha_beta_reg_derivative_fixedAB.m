% opt_method:
%   mult
%   quadratic
%   projected_grad
%   equation_solver
%   sqp
%   trust-region-reflective
%   active-set
%   interior-point
function [C isConverged] = nonlinear_alpha_beta_reg_derivative_fixedAB(I, A, B, alpha, beta, maxIterCnt, eps, alpha_reg, use_term_criteria)
    if (nargin < 9)
        use_term_criteria = true;
    end
    
    eps_nnz = 1e-12;
    minIterCnt = 50;
        
    %[A B C] = nonlinear_init_als(I, eps);
    
    C = nmf_alpha_beta_fixedA(I, A, alpha, beta, maxIterCnt, eps);
    
    %[A C] = nmf_normalize_prod(A, C);
    
    %A = rand(size(I, 1), 1);
    %C = rand(1, size(I, 2));
    %B = rand(size(I, 1), 1);
    
    %B = zeros(size(A)) + eps;
    
    reg_C_first = @(A, B, C)  2 * alpha_reg * sum(bsxfun(@times, 1 + B * C, B ./ A), 1);
    reg_C_second = @(A, B, C)  2 * alpha_reg * sum((B .^ 2) ./ A, 1);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha == 0 && beta == 0)
            F = (A * C) ./ (1 + B * C);
            C = C +  (2 * sum((1 ./ (eps_nnz + bsxfun(@plus, C, B * (C .^ 2)))) .* (log(I + eps_nnz) - log(F + eps_nnz)), 1) + reg_C_first(A, B, C)) ./ ...
                (2 * sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz)) + 1), 1) + ...
                reg_C_second(A, B, C));            
            C = max(C, eps_nnz);
        elseif (alpha == 0)
            F = (A * C) ./ (1 + B * C);
            C = C - (beta^2 * sum(power_my(F + eps_nnz, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)) .* (1 ./ (eps_nnz + A * (C .^ 2))), 1) + ...
                reg_C_first(A, B, C)) ./ ...
                (beta^2 * sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* power_my(F + eps_nnz, beta) .* ...
                (1 - (2 * B * C - beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz))), 1) + reg_C_second(A, B, C));
            C = max(C, eps_nnz);
        elseif (beta == 0)
            C = C - ((1 / alpha) * sum((power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) - power_my(eps_nnz + I, alpha)) ./ ...
                bsxfun(@times, C, 1 + B * C), 1) + reg_C_first(A, B, C)) ./ ...
                ((1 / alpha) * sum(((alpha - 1 - 2 * B * C) .* power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) + ...
                (1 + 2 * B * C) .* (power_my(I + eps_nnz, alpha))) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + reg_C_second(A, B, C));
            C = max(C, eps_nnz);
        elseif (alpha == -beta)
            direction = - 0.5 * ((1 / alpha) * sum((1 - power_my(eps_nnz + (I .* (1 + B * C)) ./ (A * C), alpha)) ./ ...
                (eps_nnz + bsxfun(@times, C, 1 + B * C)), 1) + ...
                reg_C_first(A, B, C))./ ...
                ((1 / alpha) * sum((power_my((I .* (1 + B * C)) ./ (A * C) + eps_nnz, alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + reg_C_second(A, B, C));
            C = C + direction;
            C = max(C, eps_nnz);
        else
            % optimizing C
            F = (A * C) ./ (1 + B * C);
            C = C - ((1 ./ (alpha * C .^ 2)) .* sum(bsxfun(@times, 1 ./ A, power_my(F, beta + 1) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha))), 1) + ...
                reg_C_first(A, B, C)) ./ ...
                ((1 ./ (alpha * C .^ 2)) .* sum((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* (power_my(F, alpha) .* (alpha - 2 * B * C + beta - 1) + ...
                power_my(I + eps_nnz, alpha) .* (2 * B * C - beta + 1)), 1) + ...
                reg_C_second(A, B, C));
            C = max(C, eps_nnz);
        end
             
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        
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
        %fprintf('%d: %f\n', currIter, currQuality);
        %fprintf('%d: %e\n', currIter, C(525));
        
    end
    
    isConverged = (currIter < maxIterCnt);
    
    C(isnan(C)) = 0;
end

function res = power_my(A, p)
    res = exp(p * log(A));
end