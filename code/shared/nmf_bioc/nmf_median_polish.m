function [A C isConverged] = nmf_median_polish(I, maxIterCnt, eps)

    [grand_effect, row_effect, column_effect, isConverged] = median_polish(log(I + eps), eps, maxIterCnt);
    
    C = exp(grand_effect + row_effect);
    A = exp(column_effect);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
end