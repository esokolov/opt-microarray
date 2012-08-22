r = 100; % pixels per inch
for i=1:length(alpha_range)
    alpha = alpha_range(i);
    parfor j=1:length(beta_range)
        beta = beta_range(j);
        nln_plot_probeset(inten_full_sliced{7},A_all{i,j}(inten_full_idx{7}, 7), B_all{i,j}(inten_full_idx{7}, 7), C_all{i,j}(7,:), 0);
        suplabel(['alpha=' num2str(alpha) ', beta=' num2str(beta)], 't');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 1900 1000]/r);
        print(gcf,'-dpng',sprintf('-r%d',r), ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' num2str(alpha) '_beta=' num2str(beta) '.png']);
        %    saveas(gcf, ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' int2str(alpha) '_beta=' int2str(beta) '.png'], 'png');
        saveas(gcf, ['nln_plot_probeset/probeset_' int2str(7) '_alpha=' num2str(alpha) '_beta=' num2str(beta) '.fig'], 'fig');        
    end
end