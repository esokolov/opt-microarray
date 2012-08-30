train_size = 200;
test_size = 200;

test_arrays = randperm(size(inten_full,2)-2);
test_arrays = sort(test_arrays(1:test_size));
inten_test = inten_full(:, [test_arrays end-1 end]);
inten_test_sliced = inten_full_sliced;

inten_train = zeros(size(inten_full,1),train_size+2);
inten_train_sliced = inten_full_sliced;

arrays = cell(1,50);

for i=1:50
    inten_test_sliced{i} = inten_full_sliced{i}(:,test_arrays);
    arrays{i} = pick_arrays(inten_full_sliced{i}(:,setdiff(1:(size(inten_full,2)-2), test_arrays)), train_size);    
    inten_train_sliced{i} = inten_full_sliced{i}(:,arrays{i});
    inten_train(inten_full_idx{i},:) = inten_full(inten_full_idx{i},[arrays{i} end-1 end]);
end
clear i