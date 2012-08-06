splicing_dist = zeros(1, size(inten_full_sliced, 1));
for i=1:size(inten_full_sliced, 1)
	splicing_dist(i) = mean(mean(inten_full_sliced{2}>100, 1))
end