function [Abest, Bbest, Cbest, Cbestcontrol, Abestsliced, Bbestsliced, reg_best, err_best, Avect, Bvect] = ...
    nonlinear_calibrate_frankenstein_model_exh(I_train, I_train_sliced, I_control_sliced, I_genes_inx, ...
    alpha, beta, maxIterCnt, eps)

alpha_Cs = [-Inf -15 -10 -5 -1 0 1 5 10 15];
control_errors = zeros(size(alpha_Cs));
err_eps = 1e-3;
err_best = Inf;

for i=1:length(alpha_Cs)
    fprintf('reg=%e: ', 10^alpha_Cs(i));
    
    fprintf('Factorizing I...');
    [A, B, C, Avect, Bvect, A_sliced, B_sliced] = nonlinear_calibrate_model_short(I_train, ...
        I_train_sliced, I_genes_inx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    fprintf(' Done\n');
    
    fprintf('Factorizing I_test with fixed A...');
    [C_control, control_errors(i)] = nonlinear_find_concentrations_witherror(I_control_sliced, A_sliced, B_sliced, ...
        @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    fprintf(' Done\n');
    
    if ((err_best-control_errors(i))/control_errors(i)>err_eps)
        Abest = A;
        Bbest = B;
        Cbest = C;
        Abestsliced = A_sliced;
        Bbestsliced = B_sliced;
        Cbestcontrol = C_control;
        reg_best = 10^alpha_Cs(i);
        err_best = control_errors(i);
        fprintf('1e%d: %f; %f\t%f\t%f\n', alpha_Cs(i), err_best, max(A_sliced{2}), max(B_sliced{2}), max(C(2,:)));
        %fprintf('1e%d: %f; %f\t%f\t%f\n', alpha_Cs(i), err_best, max(A_sliced{1}), max(B_sliced{1}), max(C(1,:)));
    end
    if i>1 && control_errors(i)/control_errors(i-1)>3
        fprintf('whoop error became too big');
        break
    end
end