function arrays = pick_arrays(I, arraysNum)
s = sum(I);
boundaries = quantile(s, [0 0.05 0.5 0.75 0.9 1]);
arrays = [];       
for i=1:5
    Carrays = find(s>=boundaries(i) & s<=boundaries(i+1));
    pick = randperm(length(Carrays));
    arrays = [arrays Carrays(pick(1:((i-1)*floor(length(pick)/4))))];
end

pick = randperm(length(arrays));
arrays = sort(arrays(pick(1:arraysNum)));