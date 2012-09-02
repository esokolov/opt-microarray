function [A B C C_dist_all B_dist_all A_dist_all] = nonlinear_alpha_beta_reg_derivative_dists(I, alpha, beta, maxIterCnt, eps, ...
        alpha_reg, use_term_criteria)

    eps_nnz = 1e-12;
    minIterCnt = 50;
        
    %[A B C] = nonlinear_init_als(I, eps);
    
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    
    %[A C] = nmf_normalize_prod(A, C);
    
    %A = rand(size(I, 1), 1) * 10;
    %C = rand(1, size(I, 2)) * 1000;
    %B = rand(size(I, 1), 1);
    
    B = zeros(size(A));
    
    saveIters = 50;
    
    C_dist_all = zeros(size(I,2), floor(maxIterCnt/saveIters));
    B_dist_all = zeros(size(I,1), floor(maxIterCnt/saveIters));
    A_dist_all = zeros(size(I,1), floor(maxIterCnt/saveIters));
    
    reg_B_first = @(A, B, C)  2 * alpha_reg * sum(bsxfun(@rdivide, C, A) .* (1 + B * C), 2);
    reg_B_second = @(A, B, C)  2 * alpha_reg ./ A;
    reg_C_first = @(A, B, C)  2 * alpha_reg * sum(bsxfun(@times, 1 + B * C, B ./ A), 1);
    reg_C_second = @(A, B, C)  2 * alpha_reg * sum((B .^ 2) ./ A, 1);
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha == 0 && beta == 0)
            A = exp((1 / size(I, 2)) * sum(log(I + eps_nnz) - log(bsxfun(@rdivide, C, 1 + B * C) + eps_nnz), 2));
            A = max(A, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            C = C +  (2 * sum((1 ./ (eps_nnz + bsxfun(@plus, C, B * (C .^ 2)))) .* (log(I + eps_nnz) - log(F + eps_nnz)), 1) + reg_C_first(A, B, C)) ./ ...
                (2 * sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz)) + 1), 1) + ...
                reg_C_second(A, B, C));            
            C = max(C, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            B = B -  (2 * sum(bsxfun(@rdivide, C, 1 + B * C) .* (log(I + eps_nnz) - log(F + eps_nnz)), 2) + reg_B_first(A, B, C)) ./ ...
                (2 * sum(bsxfun(@times, C .^ 2, (1 + B * C) .^ 2) .* (log(F + eps_nnz) - log(I + eps_nnz) + 1), 2) + reg_B_second(A, B, C));
            B = max(B, 0);
        elseif (alpha == 0)
            F = (A * C) ./ (1 + B * C);
            A = A - (beta^2 * (1 ./ (A + eps_nnz)) .* sum(power_my(F + eps_nnz, beta) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (beta^2 * (1 ./ (eps_nnz + A .^ 2)) .* sum(power_my(F + eps_nnz, beta) .* ((beta - 1) * log((F + eps_nnz) ./ (I + eps_nnz)) + 1), 2));
            A = max(A, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            C = C - (beta^2 * sum(power_my(F + eps_nnz, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)) .* (1 ./ (eps_nnz + A * (C .^ 2))), 1) + ...
                reg_C_first(A, B, C)) ./ ...
                (beta^2 * sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* power_my(F + eps_nnz, beta) .* ...
                (1 - (2 * B * C - beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz))), 1) + reg_C_second(A, B, C));
            C = max(C, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            B = B +  (reg_B_first(A, B, C) + ...
                beta^2 * (1 ./ (A + eps_nnz)) .* sum(power_my(F + eps_nnz, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (beta^2 * sum(bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2) .* power_my(F + eps_nnz, beta) .* ...
                ((beta + 1) * log((F + eps_nnz) ./ (I + eps_nnz)) + 1), 2) + reg_B_second(A, B, C));
            B = max(B, 0);
        elseif (beta == 0)
            F = (A * C) ./ (1 + B * C);
            A = A - ((1 ./ (alpha * A)) .* sum(power_my(F, alpha) - power_my(I, alpha), 2)) ./ ...
                ((1 ./ (alpha * A .^ 2)) .* sum((alpha - 1) * power_my(F, alpha) + power_my(I, alpha), 2));
            A = max(A, eps_nnz);
            
            C = C - ((1 / alpha) * sum((power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) - power_my(eps_nnz + I, alpha)) ./ ...
                bsxfun(@times, C, 1 + B * C), 1) + reg_C_first(A, B, C)) ./ ...
                ((1 / alpha) * sum(((alpha - 1 - 2 * B * C) .* power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) + ...
                (1 + 2 * B * C) .* (power_my(I + eps_nnz, alpha))) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + reg_C_second(A, B, C));
            C = max(C, eps_nnz);
            
            B = B -  ((1 / alpha) * sum(bsxfun(@times, C, power_my(I + eps_nnz, alpha) - power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha)) ./ ...
                (1 + B * C), 2) + reg_B_first(A, B, C)) ./ ...
                ((1 / alpha) * sum(bsxfun(@times, C .^ 2, (alpha + 1) * power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) - power_my(I + eps_nnz, alpha)) ./ ...
                ((1 + B * C) .^ 2), 2) + reg_B_second(A, B, C));
            B = max(B, 0);
        elseif (alpha == -beta)
            direction = -(size(I, 2) ./ (alpha * A) - ...
                (1 ./ (alpha * power_my(A, alpha + 1))) .* sum(power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha), 2)) ./ ...
                (- size(I, 2) ./ (alpha * (A .^ 2)) + ((alpha + 1) ./ (alpha * power_my(A, alpha + 2))) .* ...
                sum(power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha), 2));
            A = A + direction;
            A = max(A, eps_nnz);
            
            direction = - 0.5 * ((1 / alpha) * sum((1 - power_my(eps_nnz + (I .* (1 + B * C)) ./ (A * C), alpha)) ./ ...
                (eps_nnz + bsxfun(@times, C, 1 + B * C)), 1) + ...
                reg_C_first(A, B, C))./ ...
                ((1 / alpha) * sum((power_my((I .* (1 + B * C)) ./ (A * C) + eps_nnz, alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + reg_C_second(A, B, C));
            C = C + direction;
            C = max(C, eps_nnz);
            
            direction = - 0.5 * ((1 / alpha) * sum(power_my(eps_nnz + bsxfun(@rdivide, I, A), alpha) .* ...
                power_my(eps_nnz + bsxfun(@rdivide, 1 + B * C, C), alpha - 1) - ...
                bsxfun(@rdivide, C, 1 + B * C), 2) + reg_B_first(A, B, C)) ./ ...
                ((1 / alpha) * sum((alpha - 1) * (power_my(eps_nnz + bsxfun(@rdivide, I, A), alpha) .* ...
                power_my(eps_nnz + bsxfun(@rdivide, 1 + B * C, C), alpha - 2)) + ...
                bsxfun(@rdivide, C, 1 + B * C) .^ 2, 2) + reg_B_second(A, B, C));
            B = B + direction;
            B = max(B, 0);
        else
            % optimizing A
            %F = bsxfun(@rdivide, C, 1 + B * C);
            %A = (sum(((I + eps) .^ alpha) .* ((F + eps) .^ beta), 2) ./ sum((F + eps) .^ (alpha + beta), 2)) .^ (1 / alpha);
            F = (A * C) ./ (1 + B * C);
%             A = A - ((1/alpha) * (1 ./ (A + eps_nnz)) .* sum((F .^ beta) .* (F .^ alpha - (I + eps_nnz) .^ alpha), 2) - ...
%                 alpha_reg * sum(bsxfun(@rdivide, 1 + B * C, A) .^ 2, 2)) ./ ...
%                 ((1 / alpha) * ((1 ./ (A + eps_nnz)) .^ 2) .* ...
%                 sum((F .^ beta) .* ((alpha + beta - 1) * (F .^ alpha) - (beta - 1) * ((I + eps_nnz) .^ alpha)), 2) + ...
%                 2 * alpha_reg * sum(bsxfun(@rdivide, (1 + B * C) .^ 2, A .^ 3), 2));
            A = A - ((1/alpha) * (1 ./ (A + eps_nnz)) .* sum(power_my(F, beta) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha)), 2)) ./ ...
                ((1 / alpha) * ((1 ./ (A + eps_nnz)) .^ 2) .* ...
                sum(power_my(F, beta) .* ((alpha + beta - 1) * power_my(F, alpha) - (beta - 1) * power_my(I + eps_nnz, alpha)), 2));
            A = max(A, eps_nnz);

            % optimizing C
            F = (A * C) ./ (1 + B * C);
%             C = C - (sum((1 ./ (eps + A * (C .^ 2))) .* ((F + eps) .^ (beta + 1)) .* ((I + eps) .^ alpha - (F + eps) .^ alpha), 1)) ./ ...
%                (sum((1 ./ (eps + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((F + eps) .^ beta) .* (((F + eps) .^ alpha) .* (2 * B * C - alpha - beta + 1) - ...
%                ((I + eps) .^ alpha) .* (2 * B * C - beta + 1)), 1));
            C = C - ((1 ./ (alpha * C .^ 2)) .* sum(bsxfun(@times, 1 ./ A, power_my(F, beta + 1) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha))), 1) + ...
                reg_C_first(A, B, C)) ./ ...
                ((1 ./ (alpha * C .^ 2)) .* sum((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* (power_my(F, alpha) .* (alpha - 2 * B * C + beta - 1) + ...
                power_my(I + eps_nnz, alpha) .* (2 * B * C - beta + 1)), 1) + ...
                reg_C_second(A, B, C));
            C = max(C, eps_nnz);

            % optimizing B
            F = (A * C) ./ (1 + B * C);
            B = B - ((1 / alpha) * sum(bsxfun(@rdivide, C, 1 + B * C) .* power_my(F, beta) .* (power_my(I + eps_nnz, alpha) - power_my(F, alpha)), 2) + ...
                reg_B_first(A, B, C)) ./ ...
                ((1 / alpha) * sum((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* ...
                power_my(F, beta) .* ((alpha + beta + 1) * power_my(F, alpha) - (beta + 1) * power_my(I + eps_nnz, alpha)), 2) + ...
                reg_B_second(A, B, C));
%             B = B + ((1 / alpha) * (1 ./ (A + eps)) .* sum(((F + eps) .^ (beta + 1)) .* ((I + eps) .^ alpha - (F + eps) .^ alpha), 2)) ./ ...
%                ((1 / alpha) * sum(bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2) .* ((F + eps) .^ beta) .* ((beta + 1) * ((I + eps) .^ alpha) - ...
%                (alpha + beta + 1) * ((F + eps) .^ alpha)), 2));
            B = max(B, 0);
        end
        
        if mod(currIter, saveIters)==0
            [A_dist_all(:,currIter/saveIters), B_dist_all(:,currIter/saveIters), C_dist_all(:,currIter/saveIters)] = nonlinear_normalize_prod(A, B, C);
        end
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        
        if (currIter>minIterCnt && use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
            break;
        end
        prevQuality = currQuality;

        
        fprintf('%d: %f\t%e\t%e\t%e\n', currIter, currQuality,  max(C), max(B), max(A));
        %fprintf('%d: %e\n', currIter, C(912));
        
        %C_prev_iter = C;        
    end
    
    [A B C] = nonlinear_normalize_prod(A, B, C);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function res = power_my(A, p)
    res = exp(p * log(A));
end
