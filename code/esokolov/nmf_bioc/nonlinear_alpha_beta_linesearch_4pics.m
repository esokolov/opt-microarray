function [A B C isConverged qual_hist C_max_hist B_max_hist A_max_hist corr_B_hist test_qual_hist time iter_cnt] = nonlinear_alpha_beta_linesearch_4pics(I, alpha, beta, maxIterCnt, eps, ...
        alpha_A, alpha_B, alpha_C, use_term_criteria, I_test)
    if (nargin < 8)
        use_term_criteria = true;
    end
    
    eps_nnz = 1e-12;
    minIterCnt = 50;
        
    %[A B C] = nonlinear_init_als(I, eps_nnz);
    
    [A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);
    
    %[A C] = nmf_normalize_prod(A, C);
    
    %A = rand(size(I, 1), 1) * 10;
    %C = rand(1, size(I, 2)) * 1000;
    %B = rand(size(I, 1), 1);
    
    B = zeros(size(A));
    
    isConverged = 1;
    
    Aprev = A;
    Bprev = B;
    Cprev = C;
    
    %C_mean_hist = zeros(maxIterCnt, 1);
    C_max_hist = zeros(maxIterCnt, 1);
    qual_hist = zeros(maxIterCnt, 1);
    B_max_hist = zeros(maxIterCnt, 1);
    A_max_hist = zeros(maxIterCnt, 1);
    corr_B_hist = zeros(maxIterCnt, 1);
    test_qual_hist = zeros(maxIterCnt, 1);
    
    C_prev_iter = C;
    
    timer_id = tic;
    
    prevQuality = -1;
    for currIter = 1:maxIterCnt
        if (alpha == 0 && beta == 0)
            F = (A * C) ./ (1 + B * C);
            %A = exp((1 / size(I, 2)) * sum(log(I + eps_nnz) - log(bsxfun(@rdivide, C, 1 + B * C) + eps_nnz), 2));
            A = A - (alpha_A * A + (1 ./ A) .* sum(log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (alpha_A + (1 ./ (A .^ 2)) .* (sum(log((I + eps_nnz) ./ (F + eps_nnz)), 2) + size(I, 2)));
            A = max(A, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            C = C - ((1 ./ C) .* sum((1 ./ (1 + B * C)) .* log((F + eps_nnz) ./ (I + eps_nnz)), 1) + alpha_C * C) ./ ...
                ((1 ./ (C .^ 2)) .* sum((1 ./ ((1 + B * C) .^ 2)) .* ...
                (1 + (2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz))), 1) + alpha_C);
            C = max(C, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            B = B - (sum(bsxfun(@rdivide, C, 1 + B * C) .* log((I + eps_nnz) ./ (F + eps_nnz)), 2) + alpha_B * B) ./ ...
                (sum((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* (log((F + eps_nnz) ./ (I + eps_nnz)) + 1), 2) + alpha_B);
            B = max(B, 0);
        elseif (alpha == 0)
            F = (A * C) ./ (1 + B * C);
            A = A - (alpha_A * A + (1 ./ A) .* sum(power_my(F + eps_nnz, beta) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (alpha_A + (1 ./ (A .^ 2)) .* sum(power_my(F, beta) .* (1 + (beta - 1) * log((F + eps_nnz) ./ (I + eps_nnz))), 2));
            A = max(A, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            C = C - ((1 ./ (C .^ 2)) .* sum(bsxfun(@times, 1 ./ A, power_my(F, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz))), 1) + alpha_C * C) ./ ...
                ((1 ./ (C .^ 2)) .* sum((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* ...
                (1 + (beta - 1 - 2 * B * C) .* log((F + eps_nnz) ./ (I + eps_nnz))), 1) + alpha_C);
            C = max(C, eps_nnz);
            
            F = (A * C) ./ (1 + B * C);
            B = B + (alpha_B * B + (1 ./ A) .* sum(power_my(F, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)), 2)) ./ ...
                (sum((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* power_my(F + eps_nnz, beta) .* (1 + (beta + 1) * log((F + eps_nnz) ./ (I + eps_nnz))), 2) + alpha_B);
            B = max(B, 0);
        elseif (beta == 0)  % OK
            %A = (sum((I + eps_nnz) .^ alpha, 2) ./ sum((eps_nnz + bsxfun(@rdivide, C, 1 + B * C)) .^ alpha, 2)) .^ (1 / alpha);
            F = (A * C) ./ (1 + B * C);
            A = A - (alpha_A * A + (1 ./ (alpha * A)) .* sum(power_my(F, alpha) - power_my(I, alpha), 2)) ./ ...
                (alpha_A + (1 ./ (alpha * A .^ 2)) .* sum((alpha - 1) * power_my(F, alpha) + power_my(I, alpha), 2));
            A = max(A, eps_nnz);
            
            C = C - ((1 / alpha) * sum(((power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) - power_my((eps_nnz + I), alpha)) ./ ...
                bsxfun(@times, C, 1 + B * C), 1) + alpha_C * C) ./ ...
                ((1 / alpha) * sum(((alpha - 1 - 2 * B * C) .* (power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) + ...
                (1 + 2 * B * C) .* power_my((I + eps_nnz), alpha)) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + alpha_C);
            C = max(C, eps_nnz);
            
            B = B - ((1 / alpha) * sum(bsxfun(@times, C, power_my(I + eps_nnz, alpha) - power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha)) ./ ...
                (1 + B * C), 2) + alpha_B * B) ./ ...
                ((1 / alpha) * sum(bsxfun(@times, C .^ 2, (alpha + 1) * power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) - power_my(I + eps_nnz, alpha)) ./ ...
                ((1 + B * C) .^ 2), 2) + alpha_B);
            B = max(B, 0);
        elseif (alpha == -beta)     % OK
            %A = ((1 / size(I, 2)) * sum((bsxfun(@rdivide, I .* (1 + B * C), C) + eps_nnz) .^ alpha, 2)) .^ (1 / alpha);
            direction = -(alpha_A * A + size(I, 2) ./ (alpha * A) - ...
                (1 ./ (alpha * power_my(A, (alpha + 1)))) .* sum(power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha), 2)) ./ ...
                (alpha_A - size(I, 2) ./ (alpha * (A .^ 2)) + ((alpha + 1) ./ (alpha * power_my(A, (alpha + 2)))) .* ...
                sum(power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha), 2));
            %direction(A == 0 & direction < 0) = 0;
            A = A + direction;
            A = max(A, eps_nnz);
            
            direction = - 0.5 * ((1 / alpha) * sum((1 - power_my((eps_nnz + (I .* (1 + B * C)) ./ (A * C)), alpha)) ./ (eps_nnz + bsxfun(@times, C, 1 + B * C)), 1) + alpha_C * C)./ ...
                ((1 / alpha) * sum((power_my(((I .* (1 + B * C)) ./ (A * C) + eps_nnz), alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
                bsxfun(@times, C .^ 2, (1 + B * C) .^ 2), 1) + alpha_C);
            %direction(C == eps_nnz & direction < 0) = 0;
            C = C + direction;
            C = max(C, eps_nnz);
            
            direction = - 0.5 * ((1 / alpha) * sum((power_my((eps_nnz + bsxfun(@rdivide, I, A)), alpha)) .* (power_my((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)), (alpha - 1))) - ...
                bsxfun(@rdivide, C, 1 + B * C), 2) + alpha_B * B) ./ ...
                ((1 / alpha) * sum((alpha - 1) * ((power_my(eps_nnz + bsxfun(@rdivide, I, A), alpha)) .* (power_my((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)), (alpha - 2)))) + ...
                bsxfun(@rdivide, C, 1 + B * C) .^ 2, 2) + alpha_B);
            %direction(B == eps_nnz & direction < 0) = 0;
            B = B + direction;
            B = max(B, 0);
        else
            % optimizing A
            %F = bsxfun(@rdivide, C, 1 + B * C);
            %A = (sum(((I + eps_nnz) .^ alpha) .* ((F + eps_nnz) .^ beta), 2) ./ sum((F + eps_nnz) .^ (alpha + beta), 2)) .^ (1 / alpha);
            F = (A * C) ./ (1 + B * C);
            A = A - (alpha_A * A + (1/alpha) * (1 ./ (A + eps_nnz)) .* sum(power_my(F, beta) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha)), 2)) ./ ...
                (alpha_A + (1 / alpha) * ((1 ./ (A + eps_nnz)) .^ 2) .* ...
                sum(power_my(F, beta) .* ((alpha + beta - 1) * (power_my(F, alpha)) - (beta - 1) * (power_my(I + eps_nnz, alpha))), 2));
            A = max(A, eps_nnz);

            % optimizing C
            F = (A * C) ./ (1 + B * C);
            %C = C - (sum((1 ./ (eps_nnz + A * (C .^ 2))) .* ((F + eps_nnz) .^ (beta + 1)) .* ((I + eps_nnz) .^ alpha - (F + eps_nnz) .^ alpha), 1) + alpha_C * C) ./ ...
            %    (sum((1 ./ (eps_nnz + bsxfun(@times, C .^ 2, (1 + B * C) .^ 2))) .* ((F + eps_nnz) .^ beta) .* (((F + eps_nnz) .^ alpha) .* (2 * B * C - alpha - beta + 1) - ...
            %    ((I + eps_nnz) .^ alpha) .* (2 * B * C - beta + 1)), 1) + alpha_C);
            C = C - (alpha_C * C + (1 ./ (alpha * C .^ 2)) .* sum(bsxfun(@times, 1 ./ A, (power_my(F, beta + 1)) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha))), 1)) ./ ...
                (alpha_C + (1 ./ (alpha * C .^ 2)) .* sum((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* (power_my(F, alpha) .* (alpha - 2 * B * C + beta - 1) + ...
                (power_my(I + eps_nnz, alpha)) .* (2 * B * C - beta + 1)), 1));
            C = max(C, eps_nnz);

            % optimizing B
            F = (A * C) ./ (1 + B * C);
            B = B - ((1 / alpha) * sum(bsxfun(@rdivide, C, 1 + B * C) .* power_my(F, beta) .* (power_my(I + eps_nnz, alpha) - power_my(F, alpha)), 2) + alpha_B * B) ./ ...
                ((1 / alpha) * sum((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* ...
                power_my(F, beta) .* ((alpha + beta + 1) * power_my(F, alpha) - (beta + 1) * (power_my(I + eps_nnz, alpha))), 2) + alpha_B);
            %B = B + ((1 / alpha) * (1 ./ (A + eps_nnz)) .* sum(((F + eps_nnz) .^ (beta + 1)) .* ((I + eps_nnz) .^ alpha - (F + eps_nnz) .^ alpha), 2) + alpha_B * B) ./ ...
            %    ((1 / alpha) * sum(bsxfun(@rdivide, C .^ 2, (1 + B * C) .^ 2) .* ((F + eps_nnz) .^ beta) .* ((beta + 1) * ((I + eps_nnz) .^ alpha) - ...
            %    (alpha + beta + 1) * ((F + eps_nnz) .^ alpha)), 2) + alpha_B);
            B = max(B, 0);
        end
        
        %if mod(currIter, 1) == 0
        if false
            %x0 = [Aprev' Bprev' Cprev];
            %direction = [(A - Aprev)' (B - Bprev)' (C - Cprev)];
            %ls_res = fminbnd(@(x) line_func(x, x0, direction, I, alpha, beta, alpha_C, alpha_B), 0, 10, optimset('FunValCheck', 'on'));
            Abig = repmat(A, [1 size(I, 2)]);
            Bbig = repmat(B, [1 size(I, 2)]);
            Cbig = repmat(C, [size(I, 1) 1]);
            
            Fbig = repmat(A - Aprev, [1 size(I, 2)]);
            Gbig = repmat(B - Bprev, [1 size(I, 2)]);
            Hbig = repmat(C - Cprev, [size(I, 1) 1]);
            
            tau = linesearch_backtracking(Abig, Bbig, Cbig, Fbig, Gbig, Hbig, I, alpha, beta, alpha_A, alpha_B, alpha_C, 1e-4, 0.5, 1);

            %X = x0 + ls_res * direction;
            %A = X(1:size(I, 1))';
            %B = X((size(I, 1)+1):(2*size(I, 1)))';
            %C = X((2*size(I, 1) + 1):end);
            
            Anew = A + tau * (A - Aprev);
            Bnew = B + tau * (B - Bprev);
            Cnew = C + tau * (C - Cprev);

            if ((sum(Anew < 0) == 0) && (sum(Bnew < 0) == 0) && (sum(Cnew < 0) == 0))
                A = Anew;
                B = Bnew;
                C = Cnew;
            end

            Aprev = A;
            Bprev = B;
            Cprev = C;
        end
        
        % routines
        if (alpha_A == 0 && alpha_B == 0 && alpha_C == 0)
            [A, B, C] = nonlinear_normalize_prod(A, B, C);
        end
        
        Q = langmuir_func(A, B, C);
        currQuality = nmf_alpha_beta_divergence(I, Q, alpha, beta);
        currQuality_reg = nmf_alpha_beta_divergence(I, Q, alpha, beta) + 0.5 * alpha_A * sum(A .^ 2) + 0.5 * alpha_B * sum(B .^ 2) + 0.5 * alpha_C * sum(C .^ 2);
        
        %if (currIter > 1 && currQuality > prevQuality)
        %    isConverged = false;
        %    break;
        %end
        
        %if (currQuality > prevQuality && sum(C == 0) > 0)
        %    fprintf('%d\n', find(C == 0));
        %    I = I(:, setdiff(1:size(I, 2), find(C == 0)));
        %    [A B C] = nonlinear_alpha_beta(I, alpha, beta, maxIterCnt, eps_nnz, opt_method_C, opt_method_B, use_term_criteria);
        %    break;
        %end
        
        if (currIter>minIterCnt && use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
        %if (use_term_criteria && currIter > 1 && nonlinear_check_stopping_criteria_C(C, C_prev_iter, 1e-6))
            break;
        end
        prevQuality = currQuality;
        prevQuality_reg = currQuality_reg;
        %if (currIter > 10000)
        %    break;
        %end
        
        C_mean_hist(currIter) = mean(abs(C - C_prev_iter));
        C_max_hist(currIter) = max(abs(C - C_prev_iter));
        
        [A_norm, B_norm, C_norm] = nonlinear_normalize_prod(A, B, C);
        
        C_max_hist(currIter) = mean(C_norm);
        qual_hist(currIter) = currQuality;
        B_max_hist(currIter) = mean(B_norm);
        A_max_hist(currIter) = mean(A_norm);
        corr_B_hist(currIter) = corr(B_norm, quantile(I', 0.9)', 'type', 'Spearman');
        
        C_test = nonlinear_alpha_beta_fixedAB(I_test, A, B, alpha, beta, maxIterCnt, eps_nnz, alpha_C, 1);
        test_qual_hist(currIter) = nmf_alpha_beta_divergence(I_test, langmuir_func(A, B, C_test), alpha, beta);
        
        %fprintf('%d: %f\t%f\t%e\t%e\n', currIter, currQuality, currQuality_reg,  max(C), max(B));
        %fprintf('%d: %e\n', currIter, C(912));
        
        C_prev_iter = C;        
    end
    
    time = toc(timer_id);
    iter_cnt = currIter;
    
    isConverged = (currIter < maxIterCnt);
    
    %C_mean_hist = C_mean_hist(1:currIter);
    C_max_hist = C_max_hist(1:currIter);
    B_max_hist = B_max_hist(1:currIter);
    A_max_hist = A_max_hist(1:currIter);
    qual_hist = qual_hist(1:currIter);
    corr_B_hist = corr_B_hist(1:currIter);
    test_qual_hist = test_qual_hist(1:currIter);
    
    [A, B, C] = nonlinear_normalize_prod(A, B, C);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end

function f = line_func(x, x0, direction, I, alpha, beta, alpha_C, alpha_B)
    X = x0 + x * direction;
    A = X(1:size(I, 1))';
    B = X((size(I, 1)+1):(2*size(I, 1)))';
    C = X((2*size(I, 1) + 1):end);
    
    if (sum(A < 0) > 0 || sum(B < 0) > 0 || sum(C < 0) > 0)
        f = x * 1e200;
    else    
        f = nmf_alpha_beta_divergence(I, langmuir_func(A, B, C), alpha, beta) + alpha_C * sum(C .^ 2) + alpha_B * sum(B .^ 2);
    end
end

function res = power_my(A, p)
    res = exp(p * log(A));
end

function step_size = linesearch_backtracking(A, B, C, F, G, H, I, alpha, beta, alpha_A, alpha_B, alpha_C, c_param, rho, init_step)
    step_size = init_step;
    while (nmf_alpha_beta_divergence(I, ((A + step_size * F) .* (C + step_size * H)) ./ (1 + (B + step_size * G) .* (C + step_size * H)), alpha, beta) > ...
            nmf_alpha_beta_divergence(I, (A .* C) ./ (1 + B .* C), alpha, beta) + ...
            c_param * step_size * sum(sum(divergence_diff_summand_line_quad_code(step_size, A, B, C, F, G, H, I, alpha, beta, ...
            alpha_A, alpha_B, alpha_C, size(I, 1), size(I, 2)))))
        step_size = step_size * rho;
    end
end