figure;
subplot(2, 3, 1);
imagesc(alpha_range, beta_range, goodness_of_fit_fro(:, :)', [min(min(goodness_of_fit_fro(:, :))) 2e6]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('goodness of fit (frobenius)');

subplot(2, 3, 2);
imagesc(alpha_range, beta_range, goodness_of_fit_l1(:, :)', [min(min(goodness_of_fit_l1(:, :))) 1e9]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('goodness of fit (l1)');

subplot(2, 3, 3);
imagesc(alpha_range, beta_range, goodness_of_fit_weighted(:, :)', [min(min(goodness_of_fit_weighted(:, :))) 1e12]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('goodness of fit (weighted)');

subplot(2, 3, 4);
imagesc(alpha_range, beta_range, overfitting_fro(:, :)', [-1e4 1e4]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('overfitting (frobenius)');

subplot(2, 3, 5);
imagesc(alpha_range, beta_range, overfitting_l1(:, :)', [-2e5 2e5]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('overfitting (l1)');

subplot(2, 3, 6);
imagesc(alpha_range, beta_range, overfitting_weighted(:, :)', [-1e9 1e9]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
xlabel('\alpha');ylabel('\beta');
title('overfitting (weighted)');
%%
figure
subplot(2, 2, 1);
imagesc(alpha_range, beta_range, mad_A(:, :)');
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
title('median|A_1-A_2|');
xlabel('\alpha');ylabel('\beta');

subplot(2, 2, 2);
imagesc(alpha_range, beta_range, mad_B(:, :)');
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
title('median|B_1-B_2|');
xlabel('\alpha');ylabel('\beta');

% subplot(2, 2, 3);
% imagesc(alpha_range, beta_range, mad_C(:, :)', [0 500]);
% set(gca,'YDir','normal');
% colorbar('peer',gca);
% title('median|C_1-C_2|');

subplot(2, 2, 3);
imagesc(alpha_range, beta_range, concentration_dists(:, :)', [0 1000]);
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
title('mean|Cdist(C_1,C_2)|');
xlabel('\alpha');ylabel('\beta');

subplot(2, 2, 4);
imagesc(alpha_range, beta_range, corr_B(:, :)');
set(gca,'YDir','normal','XTick',[-1:5],'YTick',[-5:5]);
colorbar('peer',gca);
title('\rho(B, I_{90%})');
xlabel('\alpha');ylabel('\beta');
