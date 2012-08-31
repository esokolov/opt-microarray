function [C totalerror notConvergedCnt] = nonlinear_find_concentrations_witherror(I_sliced, A_sliced, B_sliced, factorization_func)
    sampleSize = size(I_sliced{1}, 2);
    G = length(I_sliced);
    
    notConvergedCnt = 0;
    totalerror = 0;
    
    C = zeros(G, sampleSize);
    for i = 1:G
        [C(i, :) isConverged finalerror] = factorization_func(I_sliced{i}, A_sliced{i}, B_sliced{i});
        
        notConvergedCnt = notConvergedCnt + ~isConverged;
        totalerror = totalerror + finalerror;
    end
    
    C(isnan(C)) = 0;
    
    C(isinf(C)) = max(max(C(~isinf(C))));
end
