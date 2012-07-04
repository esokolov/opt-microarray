function A = nmf_init_als_fixedC(I, C, eps)
    G = size(C, 1);
    
    A = max(eps, (I * C') * pinv(C * C'));
end