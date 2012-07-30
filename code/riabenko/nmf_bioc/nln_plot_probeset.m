function nln_plot_probeset(I,A,B,C)
figure
for i=1:size(I,1)
    subplot(5,5,i)
    scatter(C,I(i,:),'.');
    hold on
    x = 0:(max(C)/100):max(C);
    y = A(i) * x ./ (1 + B(i) * x);
    plot(x, y,'r');
    title('');
    xlim([0 max(C)]);
    ylim([0 quantile(I(i,:),0.995)]);
    box('on');
    title(['A=',num2str(A(i)),', B=',num2str(B(i))]);
end