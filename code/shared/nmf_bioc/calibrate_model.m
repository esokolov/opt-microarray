function [A, C, Avect] = calibrate_model(I, factorization_func)
    sampleSize = size(I, 2) - 2;
    probesets = unique(I(:, end));
    A = spalloc(size(I, 1), length(probesets), size(I, 1));
    Avect = zeros(size(I, 1), 1);
    C = zeros(length(probesets), sampleSize);
    for i = 1:length(probesets)
        fprintf('%d\n', i);
        idx = (I(:, end) == probesets(i));
        [aff, conc] = factorization_func(I(idx, 1:sampleSize));
        %[aff, conc] = nnmf(I(idx, 1:sampleSize), 1);
        %[aff conc] = nmf_smart(I(idx, 1:sampleSize), 1, 1000, 1e-6);
        %norm = sum(aff) / length(aff);
        norm = prod(aff) ^ (1 / length(aff));
        %norm = sum(aff) / 100;
        aff = aff / norm;
        conc = conc * norm;
        A(idx, i) = aff;
        Avect(idx) = aff;
        C(i, :) = conc;
    end
    
    Avect(isnan(Avect)) = 0;
    C(isnan(C)) = 0;
    A(isnan(A)) = 0;
    
    Avect(isinf(Avect)) = 1e10;
    C(isinf(C)) = 1e10;
    A(isinf(A)) = 1e10;
end