function res = ABdiv(P, Q, alpha, beta)
    eps = 1e-12;
    if (alpha == 0 && beta == 0)
        res = 0.5 * (log(P + eps) - log(Q + eps)) .^ 2;
        res = sum(sum(res));
    elseif (alpha == 0)
        res = (1 / beta^2) * (power_my(Q + eps, beta) .* ...
            (beta * log((Q + eps) ./ (P + eps)) - 1) + power_my(P + eps, beta));
        res = sum(sum(res));
    elseif (beta == 0)
        res = (1 / alpha^2) * (power_my(P + eps, alpha) .* ...
            (alpha * log((P + eps) ./ (Q + eps)) - 1) + power_my(Q + eps, alpha));
        res = sum(sum(res));
    elseif (alpha == -beta)
        res = (1 / alpha^2) * (alpha * log((Q + eps) ./ ...
            (P + eps)) + power_my((P + eps) ./ (Q + eps), alpha) - 1);
        res = sum(sum(res));
    else
        res = - (1 ./ (alpha * beta)) * (power_my(P + eps, alpha) .* ...
            (power_my(Q + eps, beta) - (alpha / (alpha + beta)) * power_my(P + eps, beta)) - ...
            (beta / (alpha + beta)) * power_my(Q + eps, alpha + beta));
        res = sum(sum(res));
    end
end

function res = power_my(A, p)
    res = exp(p * log(A));
end
