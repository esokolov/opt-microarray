function [A X] = nmf_init_svd(P, k, eps)
    [m n] = size(P);
    
    %A = rand(m, k);
    X = rand(k, n);
    
    A = max(eps, (P * X') * pinv(X * X'));
    X = max(eps, pinv(A' * A) * (A' * P));
    
    [A X] = nmf_normalize_l2(A, X);
end