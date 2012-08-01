figure;
subplot(5, 1, 1);
plot(A_max_hist, 'b', 'LineWidth', 2);
legend('A_{mean}', 'Location', 'NorthWest');
%title(['alpha_A = ' num2str(alpha_A) '; alpha_B = ' num2str(alpha_B) '; alpha_C = ' num2str(alpha_C)]);
grid on;
subplot(5, 1, 2);
plot(B_max_hist, 'r', 'LineWidth', 2);
legend('B_{mean}', 'Location', 'NorthWest');
grid on;
subplot(5, 1, 3);
plot(C_max_hist, 'k', 'LineWidth', 2);
legend('C_{mean}', 'Location', 'NorthWest');
grid on;
subplot(5, 1, 4);
plot(qual_hist, 'Color', [73/256 61/256 139/256], 'LineWidth', 2);
legend('Q', 'Location', 'NorthWest');
grid on;
subplot(5, 1, 5);
plot(corr_B_hist, 'Color', [46/256 139/256 87/256], 'LineWidth', 2);
legend('Corr(B, I)', 'Location', 'Best');
grid on;

set(gcf, 'Position', [0 0 1000 900]);