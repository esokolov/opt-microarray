[A C] = nmf_smart(I,1,1000,10^-15);

a_one = zeros(size(A));
for i=1:25
    [A_loo ~] = nmf_normalize(A(setdiff(1:end,i)), C);
    C_loo = nmf_smart_fixedA(I(setdiff(1:end,i),:), A_loo(setdiff(1:end,i)),1000,10^-15);
    a_one(i) = nmf_smart_fixedC(I(i,:),C_loo,1000,10^-15);
end
scatter(A,a_one)
C1 = nmf_smart_fixedA(I, a_one,1000,10^-15);
scatter(C, C1);

figure
for i=1:25
    subplot(5,5,i)
    scatter(log(I(i,:)),log(C1),'Marker','.');
    title(strcat('a=',num2str(a_one(i))))
    set(gca,'XLim',[1 12], 'YLim', [7 12])
end