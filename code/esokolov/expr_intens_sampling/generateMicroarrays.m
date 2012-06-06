for z = 1:10
    fprintf('%d\n', z);
%%
expr_x = 0:1:(2 * max(expr));
expr_y = ksdensity(expr, expr_x, 'support', 'positive');

%%
expr_cdf = pdfToCdf(expr_x, expr_y);

%%
%expr_intens_sample = zeros(10 * sum(probesetSize), 2);
%idx = 0;

%%
expr_sample = zeros(length(probesetSize), 1);
for i = 1:length(expr_sample)
    expr_sample(i) = generateFromCdf(expr_x, expr_cdf);
end

%%
[~, density, X, Y] = kde2d([expr intens]);%, (1.2 * max(max(expr), max(intens))) / 40);%, ...
                          % [0 0], [1.2*max(expr) 1.2*max(intens)]);

%%
intens_sample = cell(length(expr_sample), 1);
for geneNum = 1:length(expr_sample)
    b = find(expr_sample(geneNum) < X(1, :), 1, 'first');
    a = b - 1;
    xa = X(1, a);
    xb = X(1, b);
    ya = density(a, :);
    yb = density(b, :);
    intens_pdf = ya + ((yb - ya) ./ (xb - xa)) .* (expr_sample(geneNum) - xa);
    intens_cdf = pdfToCdf(X(1, :), intens_pdf);
    
    intens_sample{geneNum} = zeros(probesetSize(geneNum), 1);
    for probeNum = 1:probesetSize(geneNum)
        intens_sample{geneNum}(probeNum) = exp(generateFromCdf(X(1, :), intens_cdf));
        
        idx = idx + 1;
        expr_intens_sample(idx, 1) = exp(expr_sample(geneNum));
        expr_intens_sample(idx, 2) = intens_sample{geneNum}(probeNum);
    end
end

expr_sample = exp(expr_sample);
end