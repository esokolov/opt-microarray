function [A_norm B_norm C_norm] = nonlinear_normalize(A, B, C)
    D = mean(A, 1);
    %D = prod(A, 1) .^ (1 / size(A, 1));
    if (any(D == 0))
        %warning('Matrix A has zero columns');
        A_norm = A;
        B_norm = B;
        C_norm = C;
        return;
    end
    
    A_norm = A * diag(1 ./ D);
    B_norm = B * diag(1 ./ D);
    C_norm = diag(D) * C;
end