for i=1:maxLOIter-1
    if mod(i-1,5)==0
        figure
    end
    sub = mod(i-1,5)+1;
    subplot(3,5,sub)
   % hist(A_sliced{1}(:,i)-A_sliced{1}(:,i+1),27);
    stem(1:27,A_sliced{1}(:,i)-A_sliced{1}(:,i+1));
    xlabel('probe number')
    ylabel(['A_{',num2str(i),'}-A_{',num2str(i+1),'}'])
    xlim([0.5 27.5])
    ylim([-0.5 1.5])
    
    subplot(3,5,sub+5)
    %hist(B_sliced{1}(:,i)-B_sliced{1}(:,i+1),27);
    stem(1:27,B_sliced{1}(:,i)-B_sliced{1}(:,i+1));
    xlabel('probe number')
    ylabel(['B_{',num2str(i),'}-B_{',num2str(i+1),'}'])   
    xlim([0.5 27.5])
    ylim([-4e-3 1e-3]);
    
    subplot(3,5,sub+2*5)
%     x = setdiff(1:1000, find(isnan(C_sliced{1}(i+1,:)));
%     y = C_sliced{1}(i,:)-C_sliced{1}(i+1,:);
%     y = y(x);
    %hist(C_sliced{1}(i,:)-C_sliced{1}(i+1,:),100);
    stem(1:1000,C_sliced{1}(i,:)-C_sliced{1}(i+1,:));
    xlabel('array number')
    ylabel(['C_{',num2str(i),'}-C_{',num2str(i+1),'}']) 
    xlim([0.5 1000.5])
%    ylim([-4e5 1e5]);
end