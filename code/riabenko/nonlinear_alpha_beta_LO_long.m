function [A_sliced B_sliced C_sliced W_sliced isConverged] = nonlinear_alpha_beta_LO_long(inten_sliced, alpha, beta, maxIterCnt, eps, use_term_criteria)
probesetsCnt = length(inten_sliced);
arraysCnt = size(inten_sliced{1},2);
maxLOIter = 6;


W_sliced = cell(probesetsCnt, 1);
A_sliced = cell(probesetsCnt, 1);
B_sliced = cell(probesetsCnt, 1);
C_sliced = cell(probesetsCnt, 1);
isConverged = cell(probesetsCnt, 1);

for ind=1:probesetsCnt    
    W_sliced{ind} = ones(size(inten_sliced{ind},1), arraysCnt);
    A_sliced{ind} = zeros(size(inten_sliced{ind},1), maxLOIter);
    B_sliced{ind} = zeros(size(inten_sliced{ind},1), maxLOIter);
    C_sliced{ind} = zeros(maxLOIter, arraysCnt);
    isConverged{ind} = zeros(1,maxLOIter);
end

for probeset_idx = 1:probesetsCnt    
    tic;
    for LOIter=1:maxLOIter
        %fprintf('LO iteration %d; ', LOIter);        
        arrays_omit = find(mean(W_sliced{probeset_idx})<0.5);
       % fprintf('%d arrays omitted completely; ', length(arrays_omit));
        arrays_keep = setdiff(1:arraysCnt, arrays_omit);
        W_sliced{probeset_idx}(:,arrays_omit) = 0;
        if sum(sum( W_sliced{probeset_idx}(:,arrays_keep),2)==0)>0
            fprintf('EMERGENCY: probe completely omitted');
        end
        
        I = inten_sliced{probeset_idx}(:,arrays_keep);
        W = W_sliced{probeset_idx}(:,arrays_keep);
        
        [A B C isConverged{probeset_idx}(LOIter)] = nonlinear_alpha_beta_weighted(I, W, alpha, beta, maxIterCnt, eps, 0, 0, use_term_criteria);
        A_sliced{probeset_idx}(:,LOIter) = A;
        B_sliced{probeset_idx}(:,LOIter) = B;
        C_sliced{probeset_idx}(LOIter,arrays_keep) = C;
        
        if length(arrays_omit)>1
            C_sliced{probeset_idx}(LOIter,arrays_omit) = nonlinear_alpha_beta_fixedAB(inten_sliced{probeset_idx}(:,arrays_omit), A, B, alpha, beta, maxIterCnt, eps, 0, use_term_criteria);
        end       
        
        error = (langmuir_func(A,B,C)-I) ./ I .* repmat(C,size(I,1),1);
        bound = quantile(error(W==1),0.95);        
        W_sliced{probeset_idx}(:,arrays_keep) = W_sliced{probeset_idx}(:,arrays_keep).*(error<bound);                
        
       % fprintf('%d iterations, %f sec; ', step,time);
    end    
    fprintf('probeset %d: %d arrays omitted completely\n', probeset_idx, length(arrays_omit)); toc;
%    nln_plot_probeset(I,A,B,C,0)
%    suplabel([num2str(sum(sum(W_sliced{probeset_idx}))/numel(W_sliced{probeset_idx})*100), '% probes left, ', num2str(length(arrays_omit)/arraysCnt*100), '% arrays excluded completely'], 't');
%    set(gcf, 'Position', [0 0 1900 1000]);
%    saveas(gcf, ['nln_plot_probeset/probeset_' int2str(probeset_idx) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.png'], 'png');
%    saveas(gcf, ['nln_plot_probeset/probeset_' int2str(probeset_idx) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.fig'], 'fig');
end