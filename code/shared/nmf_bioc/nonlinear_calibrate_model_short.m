function [A, B, C, Avect, Bvect, A_sliced, B_sliced] = nonlinear_calibrate_model_short(I, I_sliced, I_genes_idx, ...
                                                                                    factorization_func)
sampleSize = size(I, 2) - 2;
G = length(I_sliced);

C = zeros(G, sampleSize);
A_sliced = cell(G, 1);
B_sliced = cell(G, 1);
%
%     notConvergedCnt = 0;
%     normProblemsCnt = 0;
%
%     time = zeros(G, 1);
%     iter_cnt = zeros(G, 1);


parfor i = 1:G
    %[aff, b, conc, isConverged, ~, ~, ~, ~, ~, ~, time(i), iter_cnt(i)] = factorization_func(I_sliced{i});
    [aff, b, conc] = factorization_func(I_sliced{i});
    
    %notConvergedCnt = notConvergedCnt + ~isConverged;
    
    %norm = prod(aff) ^ (1 / length(aff));
    %norm = mean(aff);
    %aff = aff / norm;
    %conc = conc * norm;
    %         if (renorm)
    %             [aff, b, conc, isProblem] = nonlinear_normalize_prod(aff, b, conc);
    %
    %             normProblemsCnt = normProblemsCnt + isProblem;
    %         end
    %
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
    
    %A_sliced{i}(isnan(A_sliced{i})) = 0;
    %A_sliced{i}(isinf(A_sliced{i})) = 1e10;
    %fprintf('%d\n', i);
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