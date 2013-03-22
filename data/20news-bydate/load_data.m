% test = importdata('test.data');
% test_label = importdata('test.label');
% test_S = sparse(test(:,1), test(:,2), test(:,3));

train = importdata('train.data');
train_label = importdata('train.label');
train_S = sparse(train(:,1), train(:,2), train(:,3));

sample_ids = [];
for i=1:20
    sample_ids = [sample_ids; find(train_label == i,5)];
end
clear i

sample_S = train_S(sample_ids, :);
sample_S = sample_S(:, sum(sample_S)>0);


