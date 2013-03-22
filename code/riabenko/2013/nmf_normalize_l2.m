function [A_norm X_norm] = nmf_normalize_l2(A, X)
    D = mean(A, 1);
    %D = prod(A, 1) .^ (1 / size(A, 1));
    if (any(D == 0))
        %warning('Matrix A has zero columns');
        A_norm = A;
        X_norm = X;
        return;
    end
    
    A_norm = A * diag(1 ./ D);
    X_norm = diag(D) * X;
    %C_norm = C;
end