function [A, C, Avect, A_sliced, notConvergedCnt, normProblemsCnt] = calibrate_model_parallel(I, I_sliced, I_genes_idx, factorization_func)
    sampleSize = size(I, 2) - 2;
    G = length(I_sliced);

    C = zeros(G, sampleSize);
    A_sliced = cell(G, 1);
    
    notConvergedCnt = 0;
    normProblemsCnt = 0;
    
    parfor i = 1:G
        [aff, conc, isConverged] = factorization_func(I_sliced{i});
        
        notConvergedCnt = notConvergedCnt + ~isConverged;
        
        %norm = prod(aff) ^ (1 / length(aff));
        %norm = mean(aff);
        %aff = aff / norm;
        %conc = conc * norm;
        [aff, conc, isProblem] = nmf_normalize_prod(aff, conc);
        
        normProblemsCnt = normProblemsCnt + isProblem;
        
        aff(isnan(aff)) = 0;
        aff(isinf(aff)) = 1e10;
        
        A_sliced{i} = aff;
        C(i, :) = conc;
        
        A_sliced{i}(isnan(A_sliced{i})) = 0;
        A_sliced{i}(isinf(A_sliced{i})) = 1e10;
    end
    
    C(isnan(C)) = 0;    
    C(isinf(C)) = 1e10;
    
    %A = spalloc(size(I, 1), length(probesets), size(I, 1));
    Avect = zeros(size(I, 1), 1);
    x = zeros(size(I, 1), 1);
    y = zeros(size(I, 1), 1);
    
    for i = 1:G
        %idx = (I(:, end) == probesets(i));
        %A(I_genes_idx{i}, i) = A_par{i};
        Avect(I_genes_idx{i}) = A_sliced{i};
        x(I_genes_idx{i}) = I_genes_idx{i};
        y(I_genes_idx{i}) = i;
    end
    
    A = sparse(x, y, Avect);
end