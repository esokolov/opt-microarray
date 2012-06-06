function cdf = pdfToCdf(x, pdf)
    cdf = cumsum(0.5 .* (x(2:end) - x(1:end-1)) .* ...
        (pdf(2:end) + pdf(1:end-1)));
    cdf = cdf / max(cdf);
    cdf = [0 cdf];
end