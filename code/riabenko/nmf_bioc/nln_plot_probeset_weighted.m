function nln_plot_probeset_weighted(I,A,B,C,W)
figure
fig_size = ceil(sqrt(size(I, 1)));
for i=1:size(I,1)
    subplot(fig_size,fig_size,i)
    hold on    
    scatter(C(W(i,:)==1),I(i,W(i,:)==1),'.');
    scatter(C(W(i,:)==0),I(i,W(i,:)==0),'g.');
    x = 0:(max(C)/1000):max(C);
    y = A(i) * x ./ (1 + B(i) * x);
    plot(x, y,'r');
    xlim([0 max(C)]);
    ylim([0 max(I(i,~isnan(C)))]);
    box('on');
    title(['A=', num2str(A(i)), ', B=', num2str(B(i))]);
end