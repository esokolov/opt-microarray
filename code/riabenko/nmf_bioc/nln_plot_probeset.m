function nln_plot_probeset(I,A,B,C,weighted)
if weighted
    weight = I - (A * C) ./ (1 + B * C);
    weight(weight<0) = sqrt(-weight(weight<0));
    weight = weight .* I;
    maxweight = quantile(weight(:),0.95);
    midweight = quantile(weight(:),0.455);
end
figure
fig_size = ceil(sqrt(size(I, 1)));
for i=1:size(I,1)
    subplot(fig_size,fig_size,i)
    hold on
    if weighted
        i
        for j=1:length(C)
            if weight(i,j)<=maxweight
                if weight(i,j)<midweight
                    plot(C(j),I(i,j),'MarkerSize',5,'Marker','.','LineStyle','none',...
                        'Color',[0 weight(i,j)/midweight 1-weight(i,j)/midweight]);
                else
                    plot(C(j),I(i,j),'MarkerSize',5,'Marker','.','LineStyle','none',...
                        'Color',[(weight(i,j)-midweight)/maxweight 1-(weight(i,j)-midweight)/maxweight 0]);
                end
            end
        end
    else
        scatter(C,I(i,:),'.');
    end
    %Csort = sort(C);
    maxC = max(C);%Csort(round(length(C)*0.99));
    x = 0:(maxC/1000):maxC;    
    y = A(i) * x ./ (1 + B(i) * x);
    plot(x, y,'r');
    title('');
    xlim([0 maxC]);
    ylim([0 max(I(i,:))]);
    box('on');
    title(['A=', num2str(A(i)), ', B=', num2str(B(i))]);
end
%saveas(f,'645.png','png');