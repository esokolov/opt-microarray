function [C notConvergedCnt] = nonlinear_find_concentrations(I_sliced, A_sliced, B_sliced, factorization_func)
    sampleSize = size(I_sliced{1}, 2);
    G = length(I_sliced);
    
    notConvergedCnt = 0;
    
    C = zeros(G, sampleSize);
    parfor i = 1:G
        [C(i, :) isConverged] = factorization_func(I_sliced{i}, A_sliced{i}, B_sliced{i});
        
        notConvergedCnt = notConvergedCnt + ~isConverged;
    end
    
    C(isnan(C)) = 0;
    
    C(isinf(C)) = max(max(C(~isinf(C))));
end
