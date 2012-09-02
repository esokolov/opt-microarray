function res = nmf_alpha_beta_divergence(I, AC, alpha, beta)
    eps = 1e-12;
    if (alpha == 0 && beta == 0)
        res = 0.5 * (log(I + eps) - log(AC + eps)) .^ 2;
        res = sum(sum(res));
    elseif (alpha == 0)
        res = (1 / beta^2) * (power_my(AC + eps, beta) .* ...
            (beta * log((AC + eps) ./ (I + eps)) - 1) + power_my(I + eps, beta));
        res = sum(sum(res));
    elseif (beta == 0)
        res = (1 / alpha^2) * (power_my(I + eps, alpha) .* ...
            (alpha * log((I + eps) ./ (AC + eps)) - 1) + power_my(AC + eps, alpha));
        res = sum(sum(res));
    elseif (alpha == -beta)
        res = (1 / alpha^2) * (alpha * log((AC + eps) ./ ...
            (I + eps)) + power_my((I + eps) ./ (AC + eps), alpha) - 1);
        res = sum(sum(res));
    else
        res = - (1 ./ (alpha * beta)) * (power_my(I + eps, alpha) .* ...
            (power_my(AC + eps, beta) - (alpha / (alpha + beta)) * power_my(I + eps, beta)) - ...
            (beta / (alpha + beta)) * power_my(AC + eps, alpha + beta));
        res = sum(sum(res));
    end
%     if (alpha == 0 && beta == 0)
%         res = 0.5 * bsxfun(@power, log(I + eps) - log(AC + eps), 2);
%         res = sum(sum(res));
%     elseif (alpha == 0)
%         res = (1 / beta^2) * (bsxfun(@times, bsxfun(@power, AC + eps, beta), ...
%             log(bsxfun(@rdivide, bsxfun(@power, AC + eps, beta) + eps, bsxfun(@power, I + eps, beta) + eps))) - ...
%             bsxfun(@power, AC + eps, beta) +...
%             bsxfun(@power, I + eps, beta));
%         res = sum(sum(res));
%     elseif (beta == 0)
%         res = (1 / alpha^2) * (bsxfun(@times, bsxfun(@power, I + eps, alpha), ...
%             log(bsxfun(@rdivide, bsxfun(@power, I + eps, alpha) + eps, bsxfun(@power, AC + eps, alpha) + eps))) - ...
%             bsxfun(@power, I + eps, alpha) +...
%             bsxfun(@power, AC + eps, alpha));
%         res = sum(sum(res));
%     elseif (alpha == -beta)
%         res = (1 / alpha^2) * ...
%             (log(bsxfun(@rdivide, bsxfun(@power, AC + eps, alpha) + eps, bsxfun(@power, I + eps, alpha) + eps)) + ...
%             bsxfun(@power, bsxfun(@rdivide, bsxfun(@power, AC + eps, alpha) + eps, bsxfun(@power, I + eps, alpha) + eps) + eps, -1) - 1);
%         res = sum(sum(res));
%     else
%         res = - (1 / (alpha * beta)) * ...
%             (bsxfun(@times, bsxfun(@power, I + eps, alpha), bsxfun(@power, AC + eps, beta)) - ...
%             (alpha / (alpha + beta)) * bsxfun(@power, I + eps, alpha + beta) - ...
%             (beta / (alpha + beta)) * bsxfun(@power, AC + eps, alpha + beta));
%         res = sum(sum(res));
%     end
end

function res = power_my(A, p)
    res = exp(p * log(A));
end
