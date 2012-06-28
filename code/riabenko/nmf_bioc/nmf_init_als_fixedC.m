function A = nmf_init_als_fixedC(I, C, eps)    
    A = max(eps, (I *C')*pinv(C * C'));
end