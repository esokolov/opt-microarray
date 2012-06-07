%% фон --- минимальная интенсивность на чипе
bg = min(I_old(:, 4:end), [], 1);

%% вычитание фона и нормализация
I = I_old;
I(:, 4:end) = max(I(:, 4:end) - repmat(bg, [size(I, 1) 1]), 0);
I(:, 4:end) = I(:, 4:end) .* repmat(normalization, [size(I, 1) 1]);

%% калибровка
probesets = unique(I(:, 3));
A = spalloc(size(I, 1), length(probesets), size(I, 1));
Avect = zeros(size(I, 1), 1);
C = zeros(length(probesets), 33);
for i = 1:length(probesets)
    fprintf('%d\n', i);
    idx = (I(:, 3) == probesets(i));
    [aff, conc] = nnmf(I(idx, 4:end), 1);
    %norm = sum(aff) / length(aff);
    %norm = prod(aff) ^ (1 / length(aff));
    norm = sum(aff);
    aff = aff / norm;
    conc = conc * norm;
    A(idx, i) = aff;
    Avect(idx) = aff;
    C(i, :) = conc;
end

%%
R = C' \ (7000000 * ones(33, 1));
i = 1:size(C, 1);
j = 1:size(C, 1);
R = sparse(i, j, R);
sum(R * C(:, 1:30), 1)

%%
R = linprog(ones(1, size(C, 1)), -speye(size(C, 1)), zeros(size(C, 1), 1), C(:, 1:30)', ...
    (7000000 * ones(30, 1)), [], [], [], optimset('MaxIter', 1000));%, 0, []);


%%
quantiles = [0 0.01 0.02 0.03 0.04 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.4 0.5];
residuals = quantiles;
for iter = 1:length(quantiles)
    q = quantiles(iter);
    fprintf('%d %f\n', iter, q);
    bg = zeros(1, 33);
    for i = 1:33
        bg(i) = quantile(I_old(3 + i, :), q);
    end
    I = I_old;
    I(:, 4:end) = max(I(:, 4:end) - repmat(bg, [size(I, 1) 1]), 0);
    I(:, 4:end) = I(:, 4:end) .* repmat(normalization, [size(I, 1) 1]);
    probesets = unique(I(:, 3));
    A = spalloc(size(I, 1), length(probesets), size(I, 1));
    %Avect = zeros(size(I, 1), 1);
    C = zeros(length(probesets), 33);
    for i = 1:length(probesets)
        if (mod(i, 1000) == 0)
            fprintf('\t%d\n', i);
        end
        idx = (I(:, 3) == probesets(i));
        [aff, conc] = nnmf(I(idx, 4:end), 1), 'options', statset('MaxIter', 1000));
        %norm = sum(aff) / length(aff);
        norm = prod(aff) ^ (1 / length(aff));
        aff = aff / norm;
        conc = conc * norm;
        A(idx, i) = aff;
        %Avect(idx) = aff;
        C(i, :) = conc;
    end
    residuals(iter) = sum(sum((I(:, 4:end) - A * C) .^ 2));
    fprintf('%f\n', residuals(iter));
end

%%
plot(quantiles, residuals, 'LineWidth', 2);
hold on;
plot(quantiles, min_res * ones(1, length(quantiles)), '--');