f = figure;
for i=1:100
    ksdensity(splicing_dist(:,i), 'support', [0-1e-12 1+1e-12]);
    xlim([0 1]);
    saveas(f, ['probe_selection\bright_dists\density',num2str(i),'.png'],'png');
end