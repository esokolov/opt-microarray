rng = -2:0.01:2;   
a = zeros(length(rng), length(rng));
for i = 1:length(rng)
    for j = 1:length(rng)
        a(i, j) = nmf_alpha_beta_divergence(1e-6, 10, rng(i), rng(j));
    end
end
imagesc(rng, rng, a, [1 1e10]);