function [Abest, Bbest, Cbest, Cbestcontrol, Abestsliced, Bbestsliced, reg_best, err_best] = ...
        nonlinear_calibrate_frankenstein_model(I_train, I_train_sliced, I_control_sliced, I_genes_inx, ...
                                                alpha, beta, maxIterCnt, eps)

alpha_Cs = [-Inf 0];%[0 10.^[-5:20]];

%control_size = size(I_control,2)-2;

control_errors = zeros(size(alpha_Cs));
%errors_type = 1;

err_eps = 1e-3;
err_best = Inf;
reg_best = 0;
i=1;
goup = true;
goon = true;
bigstep = 5;

while goon
    fprintf('reg=%e: ', 10^alpha_Cs(i));
    
    fprintf('Factorizing I...');
    %     [A, B, C, ~, ~, A_sliced, B_sliced] = nonlinear_calibrate_model_short(I_train, ...
    %         I_train_sliced, I_genes_inx, ...
    %         @(I) nonlinear_alpha_beta_linesearch(I, alpha, beta, maxIterCnt, eps, 0, 0, 10^alpha_Cs(i), 1));
    [A, B, C, ~, ~, A_sliced, B_sliced] = nonlinear_calibrate_model_short(I_train, ...
        I_train_sliced, I_genes_inx, ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    fprintf(' Done\n');
    
    fprintf('Factorizing I_test with fixed A...');
    [C_control, control_errors(i)] = nonlinear_find_concentrations_witherror(I_control_sliced, A_sliced, B_sliced, ...
        @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_LO_reg_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    %     C_control = nonlinear_find_concentrations(I_control_sliced, A_sliced, B_sliced, ...
    %         @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, eps, 0, 1));%
    fprintf(' Done\n');
    

    %nln_plot_probeset(inten_train_sliced{2},A_sliced{2},B_sliced{2},C(2,:),0)
    %suplabel(['C train, alpha_C=' num2str(10^alpha_Cs(i))], 't');
    
    %nln_plot_probeset(inten_control_sliced{2},A_sliced{2},B_sliced{2},C_control(2,:),0)
    %suplabel(['C control, alpha_C=' num2str(10^alpha_Cs(i))], 't');
%     
%     R = I_control(:,1:end-2) - langmuir_func(A, B, C_control);
%     switch errors_type
%         case 1 
%             control_errors(i) = sum(sum(R.^2)) /  control_size;
%         case 2 
%             control_errors(i)  = sum(sum(abs(R))) / control_size;
%         case 3 
%             control_errors(i) = sum(sum(loss_asymmetric(R) .* I_control(:,1:end-2) )) / control_size;
%     end        
    
    if ((err_best-control_errors(i))/control_errors(i)>err_eps)||...
            (err_best-control_errors(i)>0 && reg_best>=10^alpha_Cs(i))
        Abest = A;
        Bbest = B;
        Cbest = C;
        Abestsliced = A_sliced;
        Bbestsliced = B_sliced;
        Cbestcontrol = C_control;
        reg_best = 10^alpha_Cs(i);
        err_best = control_errors(i);
        %fprintf('1e%d: %f; %f\t%f\t%f\n', alpha_Cs(i), err_best, max(A_sliced{2}), max(B_sliced{2}), max(C(2,:)));
        fprintf('1e%d: %f; %f\t%f\t%f\n', alpha_Cs(i), err_best, max(A_sliced{1}), max(B_sliced{1}), max(C(1,:)));
    end
    
    if i>1
        if (control_errors(i)-control_errors(i-1))/control_errors(i-1)>err_eps
            goup = ~goup;
        end
        if ~goup && (abs(control_errors(i)-control_errors(1))/control_errors(1)<err_eps)
            goup = true;
        end
        if goup
            prevUps = alpha_Cs(alpha_Cs>alpha_Cs(i));
            if isempty(prevUps)
                alpha_Cs(i+1) = alpha_Cs(i) + bigstep;
            else
                if (min(prevUps)==alpha_Cs(i)+1)
                    ind = find(alpha_Cs == min(prevUps));
                    if (control_errors(ind)-control_errors(i))/control_errors(ind)>err_eps%(alpha_Cs(ind-1)<alpha_Cs(ind))
                        goon = false;
                    else
                        candidates = setdiff(alpha_Cs(i)+1:max(alpha_Cs(2:end)), prevUps);
                        if isempty(candidates)
                            goon = false;
                        else
                            alpha_Cs(i+1) = min(candidates);
                        end
                    end
                else
                    alpha_Cs(i+1) = ceil((alpha_Cs(i)+min(prevUps))/2);
                end
            end
        else
%             if abs(control_errors(i)-control_errors(1))/control_errors(1)<err_eps
%                 goon = false;
%             else
                prevDowns = alpha_Cs(alpha_Cs<alpha_Cs(i));
                prevDowns = prevDowns(~isinf(prevDowns));
                if isempty(prevDowns)                    
                    alpha_Cs(i+1) = alpha_Cs(i) - bigstep;
                else
                    if (max(prevDowns)==alpha_Cs(i)-1)
                        ind = find(alpha_Cs == max(prevDowns));
                        if (control_errors(ind)-control_errors(i))/control_errors(ind)>err_eps%(alpha_Cs(ind-1)>alpha_Cs(ind))
                            goon = false;
                        else                            
                            candidates = setdiff(min(alpha_Cs(2:end)):alpha_Cs(i)-1, prevDowns);
                            if isempty(candidates)
                                goon = false;
                            else
                                alpha_Cs(i+1) = max(candidates);
                            end
                        end
                    else
                        alpha_Cs(i+1) = floor((alpha_Cs(i)+max(prevDowns))/2);
                    end
                end
%            end
        end
    end    
    i=i+1;
end
