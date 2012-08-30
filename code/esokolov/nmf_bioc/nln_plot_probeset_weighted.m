%function nln_plot_probeset_weighted(I,A,B,C,W)
figure
fig_size = ceil(sqrt(size(I, 1)));
for i=1:size(I,1)
    subplot(fig_size,fig_size,i)
    hold on    
    scatter(C_long(arrays_omit), inten(i,arrays_omit), 'r.');
    scatter(C(W(i,arrays_keep)==1),I(i,W(i,arrays_keep)==1),'b.');
    scatter(C(W(i,arrays_keep)==0),I(i,W(i,arrays_keep)==0),'g.');    
    x = 0:(max(C)/1000):max(C);
    y = A(i) * x ./ (1 + B(i) * x);
    plot(x, y,'r');
    xlim([0 max(C)]);
    ylim([0 max(I(i,~isnan(C)))]);
    box('on');
    title(['A=', num2str(A(i)), ', B=', num2str(B(i))]);
end