function [A_long B_long C_long isConverged] = nonlinear_alpha_beta_LO(inten, alpha, beta, maxIterCnt, eps, use_term_criteria)
arraysCnt = size(inten,2);
probesCnt = size(inten,1);
maxLOIter = 6;

W = ones(size(inten));
isConverged = 0;
A_long = zeros(probesCnt,1);
B_long = zeros(probesCnt,1);
C_long = zeros(1,arraysCnt);

%tic;
for LOIter=1:maxLOIter
    %fprintf('LO iteration %d; ', LOIter);
    arrays_omit = find(mean(W)<0.5);
    % fprintf('%d arrays omitted completely; ', length(arrays_omit));
    arrays_keep = setdiff(1:arraysCnt, arrays_omit);
    W(:,arrays_omit) = 0;
    
    probes_omit = find(mean(W,2)<0.05);
    probes_keep = setdiff(1:size(inten,1), probes_omit);
    W(probes_omit,:) = 0;
        
    I = inten(probes_keep,arrays_keep);
%    W = W(:,arrays_keep);
    
    [A B C isConverged time step] = nonlinear_alpha_beta_weighted(I, W(probes_keep,arrays_keep), alpha, beta, maxIterCnt, eps, 0, 0, use_term_criteria);
    
    fprintf('%d; ',step);
    
    C_long(arrays_keep) = C;
    A_long(probes_keep) = A;
    B_long(probes_keep) = B;
    
    if ~isempty(arrays_omit)
        C_long(arrays_omit) = nonlinear_alpha_beta_fixedAB(inten(probes_keep,arrays_omit), A, B, alpha, beta, maxIterCnt, eps, 0, use_term_criteria);
    end
    
    if ~isempty(probes_omit)
        A_long(probes_omit) = 0;
        B_long(probes_omit) = 0;
    end    
     
    error = (langmuir_func(A,B,C)-I) ./ I .* repmat(C,size(I,1),1);
    bound = quantile(error(W(probes_keep,arrays_keep)==1),0.95);
    W(probes_keep,arrays_keep) = W(probes_keep,arrays_keep).*(error<=bound);
    
    %nln_plot_probeset_weighted
    % fprintf('%d iterations, %f sec; ', step,time);
end
%fprintf('probeset %d: %d arrays omitted completely\n', probeset_idx, length(arrays_omit));
%toc;
%nln_plot_probeset(I,A,B,C,0)
%suplabel([num2str(sum(sum(W_sliced{probeset_idx}))/numel(W_sliced{probeset_idx})*100), '% probes left, ', num2str(length(arrays_omit)/arraysCnt*100), '% arrays excluded completely'], 't');
%set(gcf, 'Position', [0 0 1900 1000]);
%     saveas(gcf, ['nln_plot_probeset/probeset_' int2str(probeset_idx) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.png'], 'png');
%     saveas(gcf, ['nln_plot_probeset/probeset_' int2str(probeset_idx) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.fig'], 'fig');
