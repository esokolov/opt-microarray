bg = min(inten_full(:, 1:end-2), [], 1);
inten_full(:, 1:end-2) = max(inten_full(:, 1:end-2) - repmat(bg, [size(inten_full, 1) 1]), 0);
normalization = 80 ./ median(inten_full(:, 1:end-2), 1);
inten_full(:, 1:end-2) = inten_full(:, 1:end-2) .* repmat(normalization, [size(inten_full, 1) 1]);