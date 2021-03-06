function [C_long isConverged finalerror] = nonlinear_alpha_beta_LO_reg_fixedAB(inten, A, B, alpha, beta, maxIterCnt, eps, alpha_C, use_term_criteria)
arraysCnt = size(inten,2);
probesCnt = size(inten,1);
maxLOIter = 6;

W = ones(size(inten));
    W(inten>5000 & inten > 10 * repmat(median(inten,2), 1, arraysCnt)) = 0;
isConverged = 0;
C_long = zeros(1,arraysCnt);

%tic;
for LOIter=1:maxLOIter
    %fprintf('LO iteration %d; ', LOIter);
    arrays_omit = find(sum(W)<3);
    % fprintf('%d arrays omitted completely; ', length(arrays_omit));
    arrays_keep = setdiff(1:arraysCnt, arrays_omit);
    W(:,arrays_omit) = 0;
    
    probes_omit = find(mean(W,2)<0.05);
    probes_keep = setdiff(1:probesCnt, probes_omit);
    W(probes_omit,:) = 0;
    
    I = inten(probes_keep,arrays_keep);
    %    W = W(:,arrays_keep);
        
    [C isConverged time step] = nonlinear_alpha_beta_weighted_fixedAB(I, W(probes_keep,arrays_keep), A(probes_keep), B(probes_keep), alpha, beta, maxIterCnt, eps, alpha_C, use_term_criteria);
    
    %fprintf('%d; ',step);
    
    C_long(arrays_keep) = C;
    
    if ~isempty(arrays_omit)
        C_long(arrays_omit) = nonlinear_alpha_beta_fixedAB(inten(probes_keep,arrays_omit), A(probes_keep), B(probes_keep), alpha, beta, maxIterCnt, eps, alpha_C, use_term_criteria);
    end

    error = (langmuir_func(A(probes_keep),B(probes_keep),C)-I) ./ I .* repmat(C,size(I,1),1);
    bound = quantile(error(W(probes_keep,arrays_keep)==1),0.95);
    if bound>0
        W(probes_keep,arrays_keep) = W(probes_keep,arrays_keep).*(error<=bound);
    else
        %fprintf('ALL INTENSITIES ARE UNDERESTIMATED');
        break;
    end
    %nln_plot_probeset_weighted
    % fprintf('%d iterations, %f sec; ', step,time);
end
finalerror = sum(sum(inten .* abs(langmuir_func(A,B,C_long)-inten).*W)) / sum(sum( inten .* W));