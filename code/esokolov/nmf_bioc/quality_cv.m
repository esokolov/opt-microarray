%% кросс-валидация по разбиению выборки чипов, критерий стабильности
sampleSize = size(inten, 2) - 2;
partitionSize = fix(sampleSize / 2);
partition = datasample(1:sampleSize, partitionSize, 'Replace', false);

sample1 = inten(:, [partition (size(inten, 2)-1:end)]);
sample2 = inten(:, [setdiff(1:sampleSize, partition) (size(inten, 2)-1:end)]);

[A1, C1, Avect1] = calibrate_model(sample1);
[A2, C2, Avect2] = calibrate_model(sample2);

%quality_arrays = sum((Avect1 - Avect2) .^ 2);
tmp = Avect1 - Avect2;
tmp = tmp(~isnan(tmp));
tmp = tmp(~isinf(tmp));
quality_arrays_mad = mad(tmp, 1);
quality_arrays_ouliers = sum(abs(tmp) > 3 * std(tmp));

%% кросс-валидация по разбиению проб, критерий стабильности
sampleSize = size(inten, 2) - 2;
probesets = unique(inten(:, end));
partition = zeros(size(inten, 1), 1);
partitionCompl = zeros(size(inten, 1), 1);
partitionSize = 0;
partitionComplSize = 0;
for i = 1:length(probesets)
    idx = find(inten(:, end) == probesets(i));
    probesetSize = length(idx);
    if (probesetSize == 1)
        continue;
    end
    currPartitionSize = fix(probesetSize / 2);
    
    currPartition = datasample(idx, currPartitionSize, 'Replace', false);
    partition((partitionSize + 1):(partitionSize + currPartitionSize)) = currPartition;
    partitionSize = partitionSize + currPartitionSize;
    
    currPartitionCompl = setdiff(idx, currPartition);
    partitionCompl((partitionComplSize + 1):(partitionComplSize + length(currPartitionCompl))) = ...
        currPartitionCompl;
    partitionComplSize = partitionComplSize + length(currPartitionCompl);
end
partition = partition(1:partitionSize);
partitionCompl = partitionCompl(1:partitionComplSize);

sample1 = inten(partition, :);
sample2 = inten(partitionCompl, :);

[A1, C1, Avect1] = calibrate_model(sample1);
[A2, C2, Avect2] = calibrate_model(sample2);


tmp = C1(:) - C2(:);
tmp = tmp(~isnan(tmp));
tmp = tmp(~isinf(tmp));
quality_probes_mad = mad(tmp, 1);
quality_probes_ouliers = sum(abs(tmp) > 3 * std(tmp));