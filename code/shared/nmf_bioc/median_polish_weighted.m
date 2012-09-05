function [grand_effect, row_effect, column_effect, isConverged] = median_polish_weighted(X, W, tol, maxiter)

    X(size(X, 1) + 1, size(X, 2) + 1) = 0;
    
    prevX = X;

    for it=1:maxiter    
        % Row Sweep
        X_w = X(1:end-1, 1:end-1);
        X_w(W == 0) = nan;
        X_w = [X_w; X(end, 1:end-1)];
        rm = nanmedian(X_w, 2);
        X(:, 1:end-1) = bsxfun(@minus, X(:, 1:end-1), rm);
        X(:, end) = X(:, end) + rm;

        % Column Sweep
        X_w = X(1:end-1, 1:end-1);
        X_w(W == 0) = nan;
        X_w = [X_w X(1:end-1, end)];
        cm = nanmedian(X_w, 1);
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
