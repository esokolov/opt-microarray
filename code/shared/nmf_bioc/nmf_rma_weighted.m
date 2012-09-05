function [A C isConverged] = nmf_rma_weighted(I, W)

    %C = rmasummary((0:size(I, 1)-1)', I, 'Output', 'natural');
    %A = zeros(size(I, 1), 1);
    
    I = log(I);
    
    [grand_effect, row_effect, column_effect] = median_polish_weighted(I, W, 1e-5, 1000);
    
    A = exp(column_effect);
    C = exp(grand_effect + row_effect);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
    
    [A C] = nmf_normalize_prod(A, C);
    
    isConverged = 1;
end