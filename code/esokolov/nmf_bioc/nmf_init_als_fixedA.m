function C = nmf_init_als_fixedA(I, A, eps)
    G = size(A, 2);
    
    C = max(eps, pinv(A' * A) * (A' * I));
end