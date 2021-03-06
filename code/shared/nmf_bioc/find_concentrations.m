function [C notConvergedCnt] = find_concentrations(I_sliced, A_sliced, factorization_func)
    sampleSize = size(I_sliced{1}, 2);
    G = length(I_sliced);
    
    notConvergedCnt = 0;
    
    C = zeros(G, sampleSize);
    parfor i = 1:G
        %fprintf('%d\n', i);
        %idx = (I(:, end) == probesets(i));
        
        %C(i, :) = factorization_func(I(idx, 1:sampleSize), full(A(idx, i)));
        [C(i, :) isConverged] = factorization_func(I_sliced{i}, A_sliced{i});
        
        notConvergedCnt = notConvergedCnt + ~isConverged;
    end
    
    C(isnan(C)) = 0;
    
    C(isinf(C)) = 1e10;
end
