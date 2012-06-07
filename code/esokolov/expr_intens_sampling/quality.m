%% init
alpha = [0, 0, 0, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25, ...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, ...
    0.9, 0.9, 0.9, 0.95, 0.95, 0.95, 1, 1, 1];
alpha_unique = unique(alpha);
repl_cnt = [3, 3, 3, 3, 9, 3, 3, 3, 3];
mask = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, ...
    6, 6, 6, 7, 7, 7, 8, 8, 8, 9, 9, 9];

%% variability
variability = zeros(size(C, 1), 1);
idx = 0;
for i = 1:9
    variability = variability + sum((C(:, mask == i) - ...
        repmat(mean(C(:, mask == i), 2), [1 repl_cnt(i)])) .^ 2, 2);
end
variability = variability / 27;

%% linearity
linearity = zeros(size(C, 1), 1);
for i = 1:length(linearity)
    w = regress(C(i, :)', [alpha' ones(33, 1)]);
    predicted = [alpha_unique' ones(9, 1)] * w;
    linearity(i) = sum(repl_cnt' .* ((accumarray(mask', C(i, :)', [], @mean) - ...
        predicted) .^ 2));
end

%% FC accuracy
brain_idx = ismember(probesets, brain);
fc_accuracy_brain = sum((C(brain_idx, :) ./ repmat(C(brain_idx, end), [1 33]) - ...
    repmat(1 - alpha, [sum(brain_idx) 1])) .^ 2, 2);

heart_idx = ismember(probesets, heart);
fc_accuracy_heart = sum((C(heart_idx, :) ./ repmat(C(heart_idx, 1), [1 33]) - ...
    repmat(alpha, [sum(heart_idx) 1])) .^ 2, 2);