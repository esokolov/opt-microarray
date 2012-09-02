function [grand_effect, row_effect, column_effect, isConverged] = median_polish(X, tol, maxiter)

    X(size(X, 1) + 1, size(X, 2) + 1) = 0;
    
    prevX = X;

    for it=1:maxiter    
        % Row Sweep
        rm = nanmedian(X(:, 1:end-1), 2);
        X(:, 1:end-1) = bsxfun(@minus, X(:, 1:end-1), rm);
        X(:, end) = X(:, end) + rm;

        % Column Sweep
        cm = nanmedian(X(1:end-1, :), 1);
        X(1:end-1, :) = bsxfun(@minus, X(1:end-1, :), cm);
        X(end, :) = X(end, :) + cm;

        d = abs(abs(X - prevX));
        if (d < tol)
            break;
        end
        prevX = X;
    end

    isConverged = (it < maxiter);

    column_effect = X(1:end-1, end);
    row_effect = X(end, 1:end-1);
    grand_effect = X(end, end);

%     tmp = median(row_effect);
%     row_effect = row_effect - tmp;
%     grand_effect = grand_effect + tmp;
% 
%     tmp = median(column_effect);
%     column_effect = column_effect - tmp;
%     grand_effect = grand_effect + tmp;
end
