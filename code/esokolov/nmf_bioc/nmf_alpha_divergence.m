function res = nmf_alpha_divergence(I, AC, alpha)
    if (alpha == 0)
        T = bsxfun(@times, AC, log(bsxfun(@rdivide, AC, I)));
        res = sum(sum(T - AC + I));
    elseif (alpha == 1)
        T = bsxfun(@times, I, log(bsxfun(@rdivide, I, AC)));
        res = sum(sum(T - I + AC));
    else
        T = bsxfun(@times, bsxfun(@power, I, alpha), ...
                   bsxfun(@power, AC, 1 - alpha)) - ...
            alpha * I + (alpha - 1) * AC;
        res = sum(sum(T)) / (alpha * (alpha - 1));
    end
end