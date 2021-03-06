function [A B C isConverged time step] = nonlinear_alpha_beta_weighted(I, W, alpha, beta, maxIterCnt, eps, alpha_C, alpha_B, use_term_criteria)
if (nargin < 9)
    use_term_criteria = true;
end
tic;
%[A B C] = nonlinear_init_als(I, eps);

[A C] = nmf_alpha_beta_weighted(I, W, 1, alpha, beta, maxIterCnt, eps);
%[A C] = nmf_alpha_beta(I, 1, alpha, beta, maxIterCnt, eps);

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
        F = (A * C) ./ (1 + B * C);
        %A = exp((1 / size(I, 2)) * sum(log(I + eps_nnz) - log(bsxfun(@rdivide, C, 1 + B * C) + eps_nnz), 2));
        A = A - ((1 ./ A) .* sum((log((F + eps_nnz) ./ (I + eps_nnz))).*W, 2)) ./ ...
            ((1 ./ (A .^ 2)) .* (sum((log((I + eps_nnz) ./ (F + eps_nnz))).*W, 2) + size(I, 2)));
        A = max(A, eps_nnz);
        
        F = (A * C) ./ (1 + B * C);
        C = C - ((1 ./ C) .* sum(((1 ./ (1 + B * C)) .* log((F + eps_nnz) ./ (I + eps_nnz))).*W, 1) + alpha_C * C) ./ ...
            ((1 ./ (C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* ...
            (1 + (2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz)))).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
        
        F = (A * C) ./ (1 + B * C);
        B = B - (sum((bsxfun(@rdivide, C, 1 + B * C) .* log((I + eps_nnz) ./ (F + eps_nnz))).*W, 2) + alpha_B * B) ./ ...
            (sum(((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* (log((F + eps_nnz) ./ (I + eps_nnz)) + 1)).*W, 2) + alpha_B);
        B = max(B, 0);
    elseif (alpha == 0)
        F = (A * C) ./ (1 + B * C);
        A = A - ((1 ./ A) .* sum((power_my(F + eps_nnz, beta) .* log((F + eps_nnz) ./ (I + eps_nnz))).*W, 2)) ./ ...
            ((1 ./ (A .^ 2)) .* sum((power_my(F, beta) .* (1 + (beta - 1) * log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 2));
        A = max(A, eps_nnz);
        
        F = (A * C) ./ (1 + B * C);
        C = C - ((1 ./ (C .^ 2)) .* sum((bsxfun(@times, 1 ./ A, power_my(F, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 1) + alpha_C * C) ./ ...
            ((1 ./ (C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* ...
            (1 + (beta - 1 - 2 * B * C) .* log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
        
        F = (A * C) ./ (1 + B * C);
        B = B + (alpha_B * B + (1 ./ A) .* sum((power_my(F, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz))).*W, 2)) ./ ...
            (sum(((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* power_my(F + eps_nnz, beta) .* (1 + (beta + 1) * log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 2) + alpha_B);
        B = max(B, 0);
    elseif (beta == 0)
        F = (A * C) ./ (1 + B * C);
        A = A - ((1 ./ (alpha * A)) .* sum((power_my(F, alpha) - power_my(I, alpha)).*W, 2)) ./ ...
            ((1 ./ (alpha * A .^ 2)) .* sum(((alpha - 1) * power_my(F, alpha) + power_my(I, alpha)).*W, 2));
        A = max(A, eps_nnz);
        
        C = C - ((1 / alpha) * sum((((power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) - power_my((eps_nnz + I), alpha)) ./ ...
            bsxfun(@times, C, 1 + B * C)).*W, 1) + alpha_C * C) ./ ...
            ((1 / alpha) * sum((((alpha - 1 - 2 * B * C) .* (power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) + ...
            (1 + 2 * B * C) .* power_my((I + eps_nnz), alpha)) ./ ...
            bsxfun(@times, C .^ 2, (1 + B * C) .^ 2)).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
        
        B = B - ((1 / alpha) * sum((bsxfun(@times, C, power_my(I + eps_nnz, alpha) - power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha)) ./ ...
            (1 + B * C)).*W, 2) + alpha_B * B) ./ ...
            ((1 / alpha) * sum((bsxfun(@times, C .^ 2, (alpha + 1) * power_my(eps_nnz + (A * C) ./ (1 + B * C), alpha) - power_my(I + eps_nnz, alpha)) ./ ...
            ((1 + B * C) .^ 2)).*W, 2) + alpha_B);
        B = max(B, 0);
    elseif (alpha == -beta)
        %A = ((1 / size(I, 2)) * sum((bsxfun(@rdivide, I .* (1 + B * C), C) + eps_nnz) .^ alpha, 2)) .^ (1 / alpha);
        direction = -(size(I, 2) ./ (alpha * A) - ...
            (1 ./ (alpha * power_my(A, (alpha + 1)))) .* sum((power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha)).*W, 2)) ./ ...
            ( size(I, 2) ./ (alpha * (A .^ 2)) + ((alpha + 1) ./ (alpha * power_my(A, (alpha + 2)))) .* ...
            sum((power_my(bsxfun(@rdivide, I .* (1 + B * C), C), alpha)).*W, 2));
        A = A + direction;
        A = max(A, eps_nnz);
        
        direction = - 0.5 * ((1 / alpha) * sum(((1 - power_my((eps_nnz + (I .* (1 + B * C)) ./ (A * C)), alpha)) ./ (eps_nnz + bsxfun(@times, C, 1 + B * C))).*W, 1) + alpha_C * C)./ ...
            ((1 / alpha) * sum(((power_my(((I .* (1 + B * C)) ./ (A * C) + eps_nnz), alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
            bsxfun(@times, C .^ 2, (1 + B * C) .^ 2)).*W, 1) + alpha_C);
        %direction(C == eps_nnz & direction < 0) = 0;
        C = C + direction;
        C = max(C, eps_nnz);
        
        direction = - 0.5 * ((1 / alpha) * sum(((power_my((eps_nnz + bsxfun(@rdivide, I, A)), alpha)) .* (power_my((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)), (alpha - 1))) - ...
            bsxfun(@rdivide, C, 1 + B * C)).*W, 2) + alpha_B * B) ./ ...
            ((1 / alpha) * sum(((alpha - 1) * ((power_my(eps_nnz + bsxfun(@rdivide, I, A), alpha)) .* (power_my((eps_nnz + bsxfun(@rdivide, 1 + B * C, C)), (alpha - 2)))) + ...
            bsxfun(@rdivide, C, 1 + B * C) .^ 2).*W, 2) + alpha_B);
        %direction(B == eps_nnz & direction < 0) = 0;
        B = B + direction;
        B = max(B, 0);
    else
        % optimizing A
        %         F = bsxfun(@rdivide, C, 1 + B * C);
        %         A = (sum((((I + eps_nnz) .^ alpha) .* ((F + eps_nnz) .^ beta)) .* W, 2) ./ sum(((F + eps_nnz) .^ (alpha + beta)) .* W , 2)) .^ (1 / alpha);
        %         A = max(A, eps_nnz);
        F = (A * C) ./ (1 + B * C);
        A = A - ((1/alpha) * (1 ./ (A + eps_nnz)) .* sum((power_my(F, beta) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha))).*W, 2)) ./ ...
            ((1 / alpha) * ((1 ./ (A + eps_nnz)) .^ 2) .* ...
            sum((power_my(F, beta) .* ((alpha + beta - 1) * (power_my(F, alpha)) - (beta - 1) * (power_my(I + eps_nnz, alpha)))).*W, 2));
        A = max(A, eps_nnz);
        
        % optimizing C
        F = (A * C) ./ (1 + B * C);
        C = C - (alpha_C * C + (1 ./ (alpha * C .^ 2)) .* sum((bsxfun(@times, 1 ./ A, (power_my(F, beta + 1)) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha)))).*W, 1)) ./ ...
            (alpha_C + (1 ./ (alpha * C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* (power_my(F, alpha) .* (alpha - 2 * B * C + beta - 1) + ...
            (power_my(I + eps_nnz, alpha)) .* (2 * B * C - beta + 1))).*W, 1));
        C = max(C, eps_nnz);
        
        % optimizing B
        F = (A * C) ./ (1 + B * C);
        B = B - ((1 / alpha) * sum((bsxfun(@rdivide, C, 1 + B * C) .* power_my(F, beta) .* (power_my(I + eps_nnz, alpha) - power_my(F, alpha))).*W, 2) + alpha_B * B) ./ ...
            ((1 / alpha) * sum(((bsxfun(@rdivide, C, 1 + B * C) .^ 2) .* ...
            power_my(F, beta) .* ((alpha + beta + 1) * power_my(F, alpha) - (beta + 1) * (power_my(I + eps_nnz, alpha)))).*W, 2) + alpha_B);
        B = max(B, 0);
    end
    
    Q = langmuir_func(A, B, C);
    currQuality = nmf_alpha_beta_divergence(I.*W, Q.*W, alpha, beta);
    currQuality_reg = nmf_alpha_beta_divergence(I.*W, Q.*W, alpha, beta) + 0.5 * alpha_B * sum(B .^ 2) + 0.5 * alpha_C * sum(C .^ 2);

    if (currIter>minIterCnt && use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
        break;
    end
    prevQuality = currQuality;
    %fprintf('%d: %f\n', currIter, currQuality);
    %fprintf('%d: %f %f %e %e %e\n', currIter, currQuality, currQuality_reg -currQuality, max(A), max(B), max(C));
    
end

isConverged = (currIter < maxIterCnt);

[A, B, C] = nonlinear_normalize_prod(A, B, C);

A(isnan(A)) = 0;
C(isnan(C)) = 0;
time=toc;
step=currIter;
end


function res = power_my(A, p)
res = exp(p * log(A));
end
