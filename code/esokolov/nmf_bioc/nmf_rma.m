function [A C isConverged] = nmf_rma(I)

    C = rmasummary((0:size(I, 1)-1)', I, 'Output', 'natural');
    
    A = zeros(size(I, 1), 1);
    
    A(isnan(A)) = 0;
    C(isnan(C)) = 0;
    
    isConverged = 1;
end