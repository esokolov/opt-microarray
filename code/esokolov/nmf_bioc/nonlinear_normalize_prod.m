function [A_norm B_norm C_norm isProblem] = nonlinear_normalize_prod(A, B, C)
    isProblem = 0;
    %D = mean(A, 1);
    [A, B, C] = nonlinear_normalize(A, B, C);
    
    D = prod(A, 1) .^ (1 / size(A, 1));
    if (any(D == 0))
        %warning('Matrix A has zero columns');
        isProblem = 1;
        A_norm = A;
        B_norm = B;
        C_norm = C;
        return;
    end
    
    A_norm = A * diag(1 ./ D);
    C_norm = diag(D) * C;
    B_norm = B * diag(1 ./ D);
    %C_norm = C;
end