function [C_long isConverged finalerror] = nmf_rma_LO_reg_fixedA(inten, A)
arraysCnt = size(inten,2);
probesCnt = size(inten,1);
maxLOIter = 6;

W = ones(size(inten));
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
    
    [C isConverged] = nmf_rma_fixedA_weighted(I, W(probes_keep,arrays_keep), A(probes_keep));
    
    %fprintf('%d; ',step);
    
    C_long(arrays_keep) = C;
    
    if ~isempty(arrays_omit)
        C_long(arrays_omit) = nmf_rma_fixedA(inten(probes_keep,arrays_omit), A(probes_keep));
    end

    error = (A*C-I) ./ I .* repmat(C,size(I,1),1);
    bound = quantile(error(W(probes_keep,arrays_keep)==1),0.95);
    W(probes_keep,arrays_keep) = W(probes_keep,arrays_keep).*(error<=bound);
    %nln_plot_probeset_weighted
    % fprintf('%d iterations, %f sec; ', step,time);
end
finalerror = sum(sum(inten .* abs(A*C_long-inten).*W))/arraysCnt;