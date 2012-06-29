[A C] = nmf_smart(I,1,1000,10^-15);

% figure
% for i=1:25
%     subplot(5,5,i)
%     scatter(log(I(i,:)),log(C1),'Marker','.');
%     title(strcat('a=',num2str(A1(i))))
%     set(gca,'XLim',[1 12], 'YLim', [3 8])
% end
A_start = A;
C_start = C;
[A1_start C1_start] = normalize_prod(A,C);

%A_one = zeros(size(A));
maxIterCnt = 50;
As = zeros(length(A), maxIterCnt+1);
As(:,1) = A;
for currIter = 1:maxIterCnt
    for i=1:25
        A_loo = A(setdiff(1:end,i)) / sum(abs(A(setdiff(1:end,i))));
        C_loo = nmf_smart_fixedA(I(setdiff(1:end,i),:), A_loo,100,10^-6);
        A(i) = nmf_smart_fixedC(I(i,:),C_loo,100,10^-6);
        A = A / sum(abs(A));
    end
    As(:,currIter+1) = A;
end
figure
hold on
for i=1:1+maxIterCnt
    plot(1:25, As(:,i),'Color', [0 (i-1)/(1+maxIterCnt), 1-(i-1)/(1+maxIterCnt)])
end




C = nmf_smart_fixedA(I, A,1000,10^-15);
figure('Name', num2str(maxIterCnt))
subplot(2,2,1)
scatter(A_start, A); xlabel('A\_start'); ylabel('A'); 
xlim([0 max([A A_start])]);
subplot(2,2,2)
scatter(C_start, C); xlabel('C\_start'); ylabel('C');
[A1 C1] = normalize_prod(A,C);
subplot(2,2,3)
scatter(A1_start, A1); xlabel('A1\_start'); ylabel('A1');
subplot(2,2,4)
scatter(C1_start, C1); xlabel('C1\_start'); ylabel('C1');


figure
for i=1:25
    subplot(5,5,i)
    scatter(log(I(i,:)),log(C1),'Marker','.');
    title(strcat('a=',num2str(a_one(i))))
    set(gca,'XLim',[1 12], 'YLim', [7 12])
end