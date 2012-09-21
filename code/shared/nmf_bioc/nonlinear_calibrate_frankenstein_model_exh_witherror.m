function [Abest, Bbest, Cbest, Cbestcontrol, Abestsliced, Bbestsliced, reg_best, validation_err_best, train_err_best, Abestvect, Bbestvect] = ...
    nonlinear_calibrate_frankenstein_model_exh_witherror(I_train, I_train_sliced, I_control_sliced, I_genes_inx, ...
    alpha, beta, maxIterCnt, eps)

alpha_Cs = [-Inf -15 -10 -5 -2 -1 0 1 2 5];
validation_errors = zeros(size(alpha_Cs));
err_eps = 1e-3;
validation_err_best = Inf;
train_err_best = Inf;

for i=1:length(alpha_Cs)
    fprintf('reg=%e: ', 10^alpha_Cs(i));
    
    fprintf('Factorizing I...\n');
    [A, B, C, Avect, Bvect, A_sliced, B_sliced, train_err] = nonlinear_calibrate_model_witherror(I_train, ...
        I_train_sliced, I_genes_inx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    %fprintf(' Done\n');
    
    fprintf('Factorizing I_test with fixed A...\n');
    [C_control, validation_errors(i)] = nonlinear_find_concentrations_witherror(I_control_sliced, A_sliced, B_sliced, ...
        @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    %fprintf(' Done\n');
    
    if ((validation_err_best-validation_errors(i))/validation_errors(i)>err_eps)
        Abest = A;
        Bbest = B;
        Cbest = C;
        Abestvect = Avect;
        Bbestvect = Bvect;
        Abestsliced = A_sliced;
        Bbestsliced = B_sliced;
        Cbestcontrol = C_control;
        reg_best = 10^alpha_Cs(i);
        validation_err_best = validation_errors(i);
        train_err_best = train_err;
        %fprintf('1e%d: %f; %f\t%f\t%f\n', alpha_Cs(i), validation_err_best, max(A_sliced{2}), max(B_sliced{2}), max(C(2,:)));
        fprintf('1e%d: %f - %f; %f\t%f\t%f\t%f \n', alpha_Cs(i), train_err_best, validation_err_best, max(A_sliced{2}), max(B_sliced{2}), max(C(2,:)), max(C_control(2,:)));
    end
    if i>1 && validation_errors(i)/validation_errors(i-1)>1.5
        fprintf('iteration stopped - error became too big\n');
        break
    end
end