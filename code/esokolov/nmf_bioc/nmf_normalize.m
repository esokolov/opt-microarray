function [A_norm C_norm] = nmf_normalize(A, C)
    D = sum(A, 1);
    %D = prod(A, 1) .^ (1 / size(A, 1));
    if (any(D == 0))
        warning('Matrix A has zero columns');
        A_norm = A;
        C_norm = C;
        return;
    end
    
    A_norm = A * diag(1 ./ D);
    %C_norm = diag(D) * C;
    C_norm = C;
end