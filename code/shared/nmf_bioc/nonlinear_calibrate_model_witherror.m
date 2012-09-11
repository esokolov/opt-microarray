function [A, B, C, Avect, Bvect, A_sliced, B_sliced, totalerror] = nonlinear_calibrate_model_witherror(I, I_sliced, I_genes_idx, ...
                                                                                    factorization_func)
sampleSize = size(I, 2) - 2;
G = length(I_sliced);

C = zeros(G, sampleSize);
A_sliced = cell(G, 1);
B_sliced = cell(G, 1);

totalerror = 0;

parfor i = 1:G
    [aff, b, conc, finalerror] = factorization_func(I_sliced{i});
    
    totalerror = totalerror + finalerror;

    if sum(isnan([aff;b]))>0
        aff(isnan(aff)) = 0;       
        b(isnan(aff)) = 0;
    end

    if sum(isinf([aff;b]))>0
        aff(isinf(aff)) = max(max(aff(~isinf(aff))));
        b(isinf(b)) = max(max(b(~isinf(b))));
    end
    
    
    A_sliced{i} = aff;
    B_sliced{i} = b;
    C(i, :) = conc;
    
    fprintf('probeset %d finished; ', i);
end
fprintf('\n');
C(isnan(C)) = 0;
C(isinf(C)) = max(max(C(~isinf(C))));

%A = spalloc(size(I, 1), length(probesets), size(I, 1));
Avect = zeros(size(I, 1), 1);
x_a = zeros(size(I, 1), 1);
y_a = zeros(size(I, 1), 1);

Bvect = zeros(size(I, 1), 1);
x_b = zeros(size(I, 1), 1);
y_b = zeros(size(I, 1), 1);

for i = 1:G
    %idx = (I(:, end) == probesets(i));
    %A(I_genes_idx{i}, i) = A_par{i};
    Avect(I_genes_idx{i}) = A_sliced{i};
    x_a(I_genes_idx{i}) = I_genes_idx{i};
    y_a(I_genes_idx{i}) = i;
    
    Bvect(I_genes_idx{i}) = B_sliced{i};
    x_b(I_genes_idx{i}) = I_genes_idx{i};
    y_b(I_genes_idx{i}) = i;
end

A = sparse(x_a, y_a, Avect);
B = sparse(x_b, y_b, Bvect);
end