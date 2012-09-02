function [I_sliced, I_genes_idx] = slice_I(I)
    sampleSize = size(I, 2) - 2;
    probesets = unique(I(:, end));
    
    [~, idx] = sort(I(:, end));
    I = I(idx, :);
    
    I_sliced = cell(length(probesets), 1);
    I_genes_idx = cell(length(probesets), 1);
    %A_sliced = cell(length(probesets), 1);
    for i = 1:length(probesets)
        fprintf('%d\n', i);
        idx = (I(:, end) == probesets(i));
        I_sliced{i} = I(idx, 1:sampleSize);
        I_genes_idx{i} = find(idx);
        %A_sliced{i} = full(A(idx, i));
    end
end