inten = inten_full(:, [1:500 1001:1002]);
inten_sliced = inten_full_sliced;
for i = 1:length(inten_sliced)
    inten_sliced{i} = inten_sliced{i}(:, 1:500);
end

inten_test = inten_full(:, 501:end);
inten_test_sliced = inten_full_sliced;
for i = 1:length(inten_test_sliced)
    inten_test_sliced{i} = inten_test_sliced{i}(:, 501:end);
end

maxIterCnt = 1000;

alpha = -0.5;
beta = -0.5;
fprintf('Alpha = %f, Beta = %f\n', alpha, beta);

fprintf('Factorizing I without regularization...');
[A_nr, B_nr, C_nr, Avect_nr, Bvect_nr, A_sliced_nr, B_sliced_nr, notConvergedCnt_nr] = nonlinear_calibrate_model(inten, inten_sliced, ...
    inten_full_idx, @(I) nonlinear_alpha_beta_linesearch(I, alpha, beta, maxIterCnt, 1e-6, 0, 0, 0, 1), false);
fprintf(' Done\n');

fprintf('Factorizing I_test with Q-regulariztion');
[A_quad, B_quad, C_quad, Avect_quad, Bvect_quad, A_sliced_quad, B_sliced_quad, notConvergedCnt_quad] = nonlinear_calibrate_model(inten, inten_sliced, ...
    inten_full_idx, @(I) nonlinear_alpha_beta_linesearch(I, alpha, beta, maxIterCnt, 1e-6, 0, 0.01, 1e-10, 1), false);
fprintf(' Done\n');

fprintf('Factorizing I_test with V-regulariztion');
[A_voron, B_voron, C_voron, Avect_voron, Bvect_voron, A_sliced_voron, B_sliced_voron, notConvergedCnt_voron] = nonlinear_calibrate_model(inten, inten_sliced, ...
    inten_full_idx, @(I) nonlinear_alpha_beta_reg_derivative(I, alpha, beta, maxIterCnt, 1e-6, 1e-4, 1), false);
fprintf(' Done\n');

fprintf('Factorizing I_test with fixed A (no regularization)...');
[C_test_nr, notConvergedCnt_nr_fixed] = nonlinear_find_concentrations(inten_test_sliced, A_sliced, B_sliced, ...
    @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, 1e-6, 0, 1));
fprintf(' Done\n');

fprintf('Factorizing I_test with fixed A (Q_reg)...');
[C_test_quad, notConvergedCnt_quad_fixed] = nonlinear_find_concentrations(inten_test_sliced, A_sliced, B_sliced, ...
    @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, 1e-6, 1e-10, 1));
fprintf(' Done\n');

fprintf('Factorizing I_test with fixed A (V_reg)...');
[C_test_voron, notConvergedCnt_voron_fixed] = nonlinear_find_concentrations(inten_test_sliced, A_sliced, B_sliced, ...
    @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_reg_derivative_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, 1e-6, 1e-4, 1));
fprintf(' Done\n');

qual_nr = sum(sum(huber_func(I_test - langmuir_func(A_nr, B_nr, C_test_nr))));
qual_quad = sum(sum(huber_func(I_test - langmuir_func(A_quad, B_quad, C_test_quad))));
qual_voron = sum(sum(huber_func(I_test - langmuir_func(A_voron, B_voron, C_test_voron))));