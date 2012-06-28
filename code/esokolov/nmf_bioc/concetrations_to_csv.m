for i = 50:50:1000
	load(['sample_size_agr_1000_' int2str(i) '.mat']);
	M = [probes_factors_smart{2}(:) probes_factors_smart{6}(:)];
	csvwrite(['sample_size_agr_1000_' int2str(i) '.tab'], M);
end
