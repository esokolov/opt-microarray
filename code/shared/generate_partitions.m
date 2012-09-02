function [partition_arrays, partition_probes, partition_probes_compl, ...
        partition_probes_sliced, partition_probes_sliced_compl] = ...
        generate_partitions(inten, inten_sliced, inten_genes_idx)

    sampleSize = size(inten, 2) - 2;
    partitionSize = fix(sampleSize / 2);
    %partition_arrays = datasample(1:sampleSize, partitionSize, 'Replace', false);
    partition_arrays = randsample(1:sampleSize, partitionSize, false);

    %probesets = unique(inten(:, end));
    partition_probes = zeros(size(inten, 1), 1);
    partition_probes_compl = zeros(size(inten, 1), 1);
    partition_probes_sliced = cell(length(inten_sliced), 1);
    partition_probes_sliced_compl = cell(length(inten_sliced), 1);
    partitionSize = 0;
    partitionComplSize = 0;
    for i = 1:length(inten_genes_idx)
        %idx = find(inten(:, end) == probesets(i));
        idx = inten_genes_idx{i};
        probesetSize = length(idx);
        if (probesetSize == 1)
        %if (probesetSize < 4)
            continue;
        end
        currPartitionSize = fix(probesetSize / 2);

        partition_probes_sliced{i} = randsample(length(idx), currPartitionSize, false);
        %partition_probes_sliced{i} = currPartition;
        currPartition = idx(partition_probes_sliced{i});
        
        partition_probes((partitionSize + 1):(partitionSize + currPartitionSize)) = currPartition;
        partitionSize = partitionSize + currPartitionSize;        

        partition_probes_sliced_compl{i} = setdiff(1:length(idx), partition_probes_sliced{i});
        currPartitionCompl = setdiff(idx, currPartition);
        partition_probes_compl((partitionComplSize + 1):(partitionComplSize + length(currPartitionCompl))) = ...
            currPartitionCompl;
        partitionComplSize = partitionComplSize + length(currPartitionCompl);
    end
    partition_probes = partition_probes(1:partitionSize);
    partition_probes_compl = partition_probes_compl(1:partitionComplSize);
    
    %for i = length(inten_genes_idx):-1:1
    %    if (isempty(partition_probes_sliced{i}))
    %        partition_probes_sliced(i) = [];
    %        partition_probes_sliced_compl(i) = [];
    %    end
    %end
end