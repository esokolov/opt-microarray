function C = nmf_init_als_fixedA(I, A, eps)
    G = size(A, 2);
    
    C = max(eps, pinv(min(A' * A, 1e150)) * (A' * I));
end