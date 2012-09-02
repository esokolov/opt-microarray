function [A_norm C_norm isProblem] = nmf_normalize_prod(A, C)
    isProblem = 0;
    %D = mean(A, 1);
    [A C] = nmf_normalize(A, C);
    
    D = prod(A, 1) .^ (1 / size(A, 1));
    if (any(D == 0))
        %warning('Matrix A has zero columns');
        isProblem = 1;
        A_norm = A;
        C_norm = C;
        return;
    end
    
    A_norm = A * diag(1 ./ D);
    C_norm = diag(D) * C;
    %C_norm = C;
end