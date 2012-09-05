function [C isConverged] = nmf_rma_fixedA(I, A)

    I = log(I);
    A = log(A);

    I = bsxfun(@minus, I, A);
    C = median(I, 1);
    t = median(A);
    %A = A - t;
    C = C + t;
    
    C = exp(C);
    
    C(isnan(C)) = 0;
    
    isConverged = 1;
end