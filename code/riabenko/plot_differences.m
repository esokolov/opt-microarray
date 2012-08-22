d = size(C_dist_all,2);
fig_size = ceil(sqrt(d-1));

figure
for i=1:d-1
    subplot(fig_size,fig_size,i)
    hist(C_dist_all(:,i)-C_dist_all(:,i+1),50)
end

figure
for i=1:d-1
    subplot(fig_size,fig_size,i)
    hist(A_dist_all(:,i)-A_dist_all(:,i+1),25)
end

figure
for i=1:d-1
    subplot(fig_size,fig_size,i)
    hist(B_dist_all(:,i)-B_dist_all(:,i+1),25)
end

