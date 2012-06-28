function res = nmf_alpha_divergence(I, AC, alpha)
    if (alpha == 0)
        T = bsxfun(@times, AC, log(bsxfun(@rdivide, AC + 1e-10, I + 1e-10)));   % грязный хак
        res = sum(sum(T - AC + I));
    elseif (alpha == 1)
        T = bsxfun(@times, I, log(bsxfun(@rdivide, I + 1e-10, AC + 1e-10)));    % +eps - грязный хак
        T(I == 0) = 0;
        res = sum(sum(T - I + AC));
    else
        T = bsxfun(@times, bsxfun(@power, I + 1e-10, alpha), ...
                   bsxfun(@power, AC + 1e-10, 1 - alpha)) - ...
            alpha * I + (alpha - 1) * AC;
        res = sum(sum(T)) / (alpha * (alpha - 1));
    end
end