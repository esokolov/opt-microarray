load('/home/affy/model/inten_allarrays_50probesets_bgnorm.mat')
alpha_range = 0:0.5:5;
beta_range = -4:0.5:4;
probesets = 1:25;
%arrays = 1:1000;
maxIterCnt = 1000;
eps = 1e-5;
r = 100; % pixels per inch
A_all = cell(length(alpha_range),length(beta_range));
B_all = cell(length(alpha_range),length(beta_range));
C_all = cell(length(alpha_range),length(beta_range));

tmp=[];
inten_sliced = inten_train_sliced(probesets);
inten_idx = inten_full_idx(probesets);
for i = 1:length(probesets)
    inten_sliced{i} = inten_train_sliced{probesets(i)} + 1;
    tmp = [tmp; inten_idx{i}];
end
inten = inten_full(tmp, [arrays end-1 end]) + 1;
clear i tmp

goodness_of_fit_fro = zeros(length(alpha_range), length(beta_range));
goodness_of_fit_l1 = zeros(length(alpha_range), length(beta_range));
goodness_of_fit_weighted = zeros(length(alpha_range), length(beta_range));
corr_B = zeros(length(alpha_range), length(beta_range));

for ind = 1:length(alpha_range)*length(beta_range)
    tic;
    i = mod(ind-1,length(alpha_range))+1;
    alpha = alpha_range(i);
    j = round((ind-i)/length(alpha_range)) + 1;
    beta = beta_range(j);
    fprintf('Alpha = %f, Beta = %f\n', alpha, beta);
    
    [A, B, C, Avect, Bvect, A_sliced, B_sliced] = nonlinear_calibrate_model_short(inten, inten_sliced, ...
        inten_idx, @(I) nonlinear_alpha_beta_LO(I, alpha, beta, maxIterCnt, eps, 1));
    
    corr_B(i, j) = corr(Bvect, quantile(inten(:, 1:end-2)', 0.9)', 'type', 'Spearman');
    
    R = inten(:, 1:end-2) - langmuir_func(A, B, C);
    goodness_of_fit_fro(i, j) = sum(sum(R.^2)) / train_size;
    goodness_of_fit_l1(i, j)  = sum(sum(abs(R))) / train_size;
    goodness_of_fit_weighted(i, j) = sum(sum(loss_asymmetric(R) .* inten(:, 1:end-2) )) / train_size;
    
    A_all{i,j} = A;
    B_all{i,j} = B;
    C_all{i,j} = C;
    
    fprintf('Goodness of fit (fro): %e\nGoodness of fit (l1): %e\nGoodness of fit (weighted): %e\nCorr: %f\n', ...
             goodness_of_fit_fro(i, j), goodness_of_fit_l1(i, j), goodness_of_fit_weighted(i, j),  corr_B(i, j));        
    toc;
    fprintf('\n');
    save(['alpha_beta_LO_minicv_200chosenarrays_25probesets_' num2str(alpha - min(alpha_range)) '_' num2str(beta - min(beta_range)) '.mat'],...
        'alpha_range', 'beta_range', 'goodness_of_fit_fro', 'goodness_of_fit_l1', 'goodness_of_fit_weighted', 'A_all', 'B_all', 'C_all', 'corr_B');

    nln_plot_probeset(inten_sliced{7},A(inten_idx{7},7),B(inten_idx{7},7),C(7,:),0)    
%    set(gcf, 'Position', [0 0 1900 1000]);
    suplabel(['alpha=' num2str(alpha) ', beta=' num2str(beta)], 't');

    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1900 1000]/r);
    print(gcf,'-dpng',sprintf('-r%d',r), ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' num2str(alpha) '_beta=' num2str(beta) '.png']);
%    saveas(gcf, ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.png'], 'png');
    saveas(gcf, ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' num2str(alpha) '_beta=' num2str(beta) '.fig'], 'fig');
end

clear i j ind alpha beta A B C Avect Bvect