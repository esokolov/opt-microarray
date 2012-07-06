function [A B C isConverged] = nonlinear_alpha_beta_fmincon(I, alpha, beta, maxIterCnt, eps)
    %[A B C] = nonlinear_init_als(I, eps);
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    %[A C] = nmf_normalize_prod(A, C);
    B = zeros(size(A));
    
    I_true = I;
    
    [P K] = size(I);
    I = I(:);
    A = A(:);
    B = B(:);
    C = C(:);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        % optimizing A
        A = fmincon(@(A_arg) cost_and_grad_A(I, A_arg, B, C, P, K), A, [], [], [], [], 0, [], [], ...
            optimset('Algorithm', 'trust-region-reflective', 'GradObj', 'on', 'Display', 'off'));
         
        % optimizing C
        C = fmincon(@(C_arg) cost_and_grad_C(I, A, B, C_arg, P, K), C, [], [], [], [], 0, [], [], ...
            optimset('Algorithm', 'trust-region-reflective', 'GradObj', 'on', 'Display', 'off'));
        
        % optimizing B
        B = fmincon(@(B_arg) cost_and_grad_B(I, A, B_arg, C, P, K), B, [], [], [], [], 0, [], [], ...
            optimset('Algorithm', 'trust-region-reflective', 'GradObj', 'on', 'Display', 'off'));
        
        % routines
        [A, B, C] = nonlinear_normalize_prod(A, B, C);
        
        [currQuality Q] = nonlinear_alpha_beta_divergence_special(I, A, B, C, P, K);
        if (nonlinear_check_stopping_criteria(I_true, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;
        
        fprintf('%d: %f\n', currIter, currQuality);
    end
    
    [~, A, B, C] = vectors_to_matrices(I, A, B, C, P, K);
    
    isConverged = (currIter < maxIterCnt);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function [cost grad] = cost_and_grad_A(I, A, B, C, P, K)
    cost = nonlinear_alpha_beta_divergence_special(I, A, B, C, P, K);
    grad = nonlinear_alpha_beta_grad_A_special(I, A, B, C, P, K);
end

function [cost grad] = cost_and_grad_B(I, A, B, C, P, K)
    cost = nonlinear_alpha_beta_divergence_special(I, A, B, C, P, K);
    grad = nonlinear_alpha_beta_grad_B_special(I, A, B, C, P, K);
end

function [cost grad] = cost_and_grad_C(I, A, B, C, P, K)
    cost = nonlinear_alpha_beta_divergence_special(I, A, B, C, P, K);
    grad = nonlinear_alpha_beta_grad_C_special(I, A, B, C, P, K);
end

function [div Q] = nonlinear_alpha_beta_divergence_special(I, A, B, C, P, K)
    [I A B C] = vectors_to_matrices(I, A, B, C, P, K);
    Q = langmuir_func(A, B, C);
    div = sum(sum(log(I) - log(Q) + Q ./ I - 1));
end

function grad = nonlinear_alpha_beta_grad_A_special(I, A, B, C, P, K)
    [I A B C] = vectors_to_matrices(I, A, B, C, P, K);
    grad = sum( (bsxfun(@rdivide, C, I) - bsxfun(@rdivide, 1 + B * C, A)) .* bsxfun(@rdivide, 1, 1 + B * C), 2);
end

function grad = nonlinear_alpha_beta_grad_B_special(I, A, B, C, P, K)
    [I A B C] = vectors_to_matrices(I, A, B, C, P, K);
    grad = -sum( ((A * C) ./ I - (1 + B * C)) .* bsxfun(@rdivide, C, (1 + B * C) .^ 2), 2);
end

function grad = nonlinear_alpha_beta_grad_C_special(I, A, B, C, P, K)
    [I A B C] = vectors_to_matrices(I, A, B, C, P, K);
    grad = sum( (bsxfun(@rdivide, A, I) - bsxfun(@rdivide, 1 + B * C, C)) .* bsxfun(@rdivide, 1, (1 + B * C) .^ 2), 1);
end

function [I A B C] = vectors_to_matrices(I, A, B, C, P, K)
    I = reshape(I, [P K]);
    A = reshape(A, [P 1]);
    B = reshape(B, [P 1]);
    C = reshape(C, [1 K]);
end

% function grad = nonlinear_alpha_beta_grad_B(I, A, B, C, alpha, beta)
%     F = (A * C) ./ (1 + B * C);
%     D = (A * (C .^ 2)) ./ ((1 + B * C) .^ 2);
%     grad = sum(-(((I + eps) .^ alpha) .* ((F + eps) .^ (beta - 1)) .* D) + ...
%         ((F + eps) .^ (alpha + beta - 1) .* D), 2);
% end