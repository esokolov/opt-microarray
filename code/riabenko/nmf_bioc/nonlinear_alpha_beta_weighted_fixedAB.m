function [C isConverged time step] = nonlinear_alpha_beta_weighted_fixedAB(I, W, A, B, alpha, beta, maxIterCnt, eps, alpha_C, use_term_criteria)
if (nargin < 10)
    use_term_criteria = true;
end
tic;

C = nmf_alpha_beta_fixedA(I, A, alpha, beta, maxIterCnt, eps);

minIterCnt = 50;
eps_nnz = 1e-12;

isConverged = 1;

prevQuality = -1;
for currIter = 1:maxIterCnt
    if (alpha == 0 && beta == 0)
        F = (A * C) ./ (1 + B * C);
        C = C - ((1 ./ C) .* sum(((1 ./ (1 + B * C)) .* log((F + eps_nnz) ./ (I + eps_nnz))).*W, 1) + alpha_C * C) ./ ...
            ((1 ./ (C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* ...
            (1 + (2 * B * C + 1) .* log((I + eps_nnz) ./ (F + eps_nnz)))).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
    elseif (alpha == 0)
        F = (A * C) ./ (1 + B * C);
        C = C - ((1 ./ (C .^ 2)) .* sum((bsxfun(@times, 1 ./ A, power_my(F, beta + 1) .* log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 1) + alpha_C * C) ./ ...
            ((1 ./ (C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* ...
            (1 + (beta - 1 - 2 * B * C) .* log((F + eps_nnz) ./ (I + eps_nnz)))).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
     elseif (beta == 0)
        C = C - ((1 / alpha) * sum((((power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) - power_my((eps_nnz + I), alpha)) ./ ...
            bsxfun(@times, C, 1 + B * C)).*W, 1) + alpha_C * C) ./ ...
            ((1 / alpha) * sum((((alpha - 1 - 2 * B * C) .* (power_my((eps_nnz + (A * C) ./ (1 + B * C)), alpha)) + ...
            (1 + 2 * B * C) .* power_my((I + eps_nnz), alpha)) ./ ...
            bsxfun(@times, C .^ 2, (1 + B * C) .^ 2)).*W, 1) + alpha_C);
        C = max(C, eps_nnz);
    elseif (alpha == -beta)
        direction = - 0.5 * ((1 / alpha) * sum(((1 - power_my((eps_nnz + (I .* (1 + B * C)) ./ (A * C)), alpha)) ./ (eps_nnz + bsxfun(@times, C, 1 + B * C))).*W, 1) + alpha_C * C)./ ...
            ((1 / alpha) * sum(((power_my(((I .* (1 + B * C)) ./ (A * C) + eps_nnz), alpha) .* (2 * B * C + alpha + 1) - 2 * B * C - 1) ./ ...
            bsxfun(@times, C .^ 2, (1 + B * C) .^ 2)).*W, 1) + alpha_C);
        C = C + direction;
        C = max(C, eps_nnz);
    else
        F = (A * C) ./ (1 + B * C);
        C = C - (alpha_C * C + (1 ./ (alpha * C .^ 2)) .* sum((bsxfun(@times, 1 ./ A, (power_my(F, beta + 1)) .* (power_my(F, alpha) - power_my(I + eps_nnz, alpha)))).*W, 1)) ./ ...
            (alpha_C + (1 ./ (alpha * C .^ 2)) .* sum(((1 ./ ((1 + B * C) .^ 2)) .* power_my(F, beta) .* (power_my(F, alpha) .* (alpha - 2 * B * C + beta - 1) + ...
            (power_my(I + eps_nnz, alpha)) .* (2 * B * C - beta + 1))).*W, 1));
        C = max(C, eps_nnz);
     end
    
    Q = langmuir_func(A, B, C);
    currQuality = nmf_alpha_beta_divergence(I.*W, Q.*W, alpha, beta);
    currQuality_reg = nmf_alpha_beta_divergence(I.*W, Q.*W, alpha, beta) + 0.5 * alpha_C * sum(C .^ 2);

    if (currIter>minIterCnt && use_term_criteria && nonlinear_check_stopping_criteria(I, Q, currQuality, prevQuality, eps))
        break;
    end
    prevQuality = currQuality;
    %fprintf('%d: %f\n', currIter, currQuality);
    fprintf('%d: %f %f %e %e\n', currIter, currQuality, currQuality_reg -currQuality, max(C));
    
end

isConverged = (currIter < maxIterCnt);

C(isnan(C)) = 0;
time=toc;
step=currIter;
end


function res = power_my(A, p)
res = exp(p * log(A));
end
