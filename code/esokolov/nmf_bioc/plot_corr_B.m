alpha_range = -4:0.25:4;
beta_range = -4:0.25:4;

figure
%imagesc(alpha_range, beta_range, goodness_of_fit_huber', [min(goodness_of_fit_huber(:)) quantile(goodness_of_fit_huber(:), 0.80)]);
%imagesc(alpha_range, beta_range, goodness_of_fit_huber(:, :, 2)', [min(min(goodness_of_fit_huber(:, :, 2))) 3*10^7]);
imagesc(alpha_range, beta_range, corr_B(:, :, 2)', [-1 1]);
set(gca,'YDir','normal');

%%
figure;
i = 10;
subplot(2, 2, 1);
imagesc(alpha_range, beta_range, corr_B_quad(:, :, i)', [-1 1]);
set(gca,'YDir','normal');
colorbar('peer',gca);
subplot(2, 2, 2);
imagesc(alpha_range, beta_range, goodness_of_fit_huber_quad(:, :, i)', [min(min(goodness_of_fit_huber_quad(:, :, i))) 3*10^7]);
set(gca,'YDir','normal');
colorbar('peer',gca);
subplot(2, 2, 3);
imagesc(alpha_range, beta_range, corr_B(:, :, i)', [-1 1]);
set(gca,'YDir','normal');
colorbar('peer',gca);
subplot(2, 2, 4);
imagesc(alpha_range, beta_range, goodness_of_fit_huber(:, :, i)', [min(min(goodness_of_fit_huber(:, :, i))) 3*10^7]);
set(gca,'YDir','normal');
colorbar('peer',gca);