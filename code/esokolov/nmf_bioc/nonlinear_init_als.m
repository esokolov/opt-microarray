function [A B C] = nonlinear_init_als(I, eps)
    G = 1;
    [P K] = size(I);
    
    A = rand(P, G);
    C = rand(G, K);
    
    A = max(eps, (I * C') * pinv(C * C'));
    C = max(eps, pinv(A' * A) * (A' * I));
    
    [A C] = nmf_normalize(A, C);
    
    B = zeros(P, G);
end