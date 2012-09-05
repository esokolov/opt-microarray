function qual = quality_functional(inten_sliced, A_sliced, B_sliced, C)
    qual = 0;
    sample_size = size(inten_sliced{1}, 2);

    for gene_idx = 1:length(inten_sliced)
        I = inten_sliced{gene_idx};
        error = (langmuir_func(A_sliced{gene_idx}, B_sliced{gene_idx}, C(gene_idx, :)) - I) ./ ...
            I .* repmat(C(gene_idx, :), size(I, 1), 1);
        bound = quantile(error(:), 0.75);
        W = ones(size(I));
        W = W .* (error <= bound);
        qual = qual + ...
            sum(sum(I .* abs(langmuir_func(A_sliced{gene_idx}, B_sliced{gene_idx}, C(gene_idx, :)) - I) .* W)) / sample_size;
    end
end