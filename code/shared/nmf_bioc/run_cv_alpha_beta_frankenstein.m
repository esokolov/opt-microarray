alpha_range = -2:0.5:4;
beta_range = -4:0.5:4;
%alpha_range = -2:0.25:2;
%beta_range = -2:0.25:4;

rep_a = zeros(length(alpha_range), length(beta_range));
rep_b = zeros(length(alpha_range), length(beta_range));
C_loo = zeros(length(alpha_range), length(beta_range));
train_error = zeros(length(alpha_range), length(beta_range));
validation_error = zeros(length(alpha_range), length(beta_range));
test_error = zeros(length(alpha_range), length(beta_range));
overfitting = zeros(length(alpha_range), length(beta_range));
reg_best = zeros(length(alpha_range), length(beta_range));

A_all = cell(length(alpha_range), length(beta_range));
B_all = cell(length(alpha_range), length(beta_range));
C_all = cell(length(alpha_range), length(beta_range));
A_sliced_all = cell(length(alpha_range), length(beta_range));
B_sliced_all = cell(length(alpha_range), length(beta_range));

maxIterCnt = 1000;
eps = 1e-5;

for i = 1:length(alpha_range)
    for j = 1:length(beta_range)
        alpha = alpha_range(i);
        beta = beta_range(j);
        fprintf('Alpha = %f, Beta = %f\n', alpha, beta);
        
        if ((alpha < -beta)&&alpha>=0) ||(alpha<0 && beta<0)
            continue;
        end
        
        [A, B, C, C_val, A_sliced, B_sliced, reg_best(i, j), train_error(i, j), validation_error(i, j), Avect, Bvect] = ...
            nonlinear_calibrate_frankenstein_model_exh_witherror(inten_train, inten_train_sliced, inten_validation_sliced, inten_full_idx, ...
            alpha, beta, maxIterCnt, eps);

        A_all{i, j} = A;
        B_all{i, j} = B;
        C_all{i, j} = C;
        
        A_sliced_all{i, j} = A_sliced;
        B_sliced_all{i, j} = B_sliced;                
        
        [C_test, test_error(i,j)] = nonlinear_find_concentrations_witherror(inten_test_sliced, A_sliced, B_sliced, ...
            @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, ...
            alpha, beta, maxIterCnt, eps, reg_best(i,j), 1));
                
        overfitting(i, j) = test_error(i, j) - train_error(i, j);
        
        C_loo(i, j) = nonlinear_get_LOO_error(inten_train_sliced, A_sliced, B_sliced, C, reg_best(i, j), ...
            alpha, beta, maxIterCnt, eps);        
        
        [~, ~, ~, Avect2, Bvect2] = nonlinear_calibrate_model_witherror(inten_test, inten_test_sliced, inten_full_idx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, reg_best(i,j), 1));
        
        rep_a(i, j) = mean((abs(Avect - Avect2)./(Avect + Avect2))>0.5);
        
        rep_b(i, j) = mean((abs(Bvect - Bvect2)./(Bvect + Bvect2))>0.5);
        
        save('cv_frankenstein_res.mat', 'rep_a', 'rep_b', 'C_loo', 'reg_best',...
            'train_error', 'validation_error', 'test_error', 'overfitting', 'A_all', 'B_all', 'C_all', ...
            'A_sliced_all', 'B_sliced_all', 'reg_best');
    end
end