alpha = 4;
beta = -3;
maxIterCnt=1000;
eps = 1e-5;

alpha_Cs = [-Inf 0];%[0 10.^[-5:20]];

test_errors_fro = zeros(size(alpha_Cs));
%test_errors_l1 = zeros(size(alpha_Cs));
%test_errors_weighted = zeros(size(alpha_Cs));

err_best = Inf;
reg_best = 0;
i=1;
goup = true;
goon = true;
while goon
    fprintf('reg=%e: ', 10^alpha_Cs(i));
    [A, B, C, Avect, Bvect, A_sliced, B_sliced] = nonlinear_calibrate_model_short(inten_train(1:32,:), ...
        inten_train_sliced(1:2), inten_full_idx(1:2), ...
        @(I) nonlinear_alpha_beta_LO_reg(I, alpha, beta, maxIterCnt, eps, 10^alpha_Cs(i), 1));
    C_test = nonlinear_find_concentrations(inten_test_sliced(1:2), A_sliced, B_sliced, ...
        @(I_arg, A_arg, B_arg) nonlinear_alpha_beta_fixedAB(I_arg, A_arg, B_arg, alpha, beta, maxIterCnt, eps, 0, 1));
    
    %nln_plot_probeset(inten_train_sliced{2},A_sliced{2},B_sliced{2},C(2,:),0)
    %suplabel(['C train, alpha_C=' num2str(10^alpha_Cs(i))], 't');
    
    %nln_plot_probeset(inten_test_sliced{2},A_sliced{2},B_sliced{2},C_test(2,:),0)
    %suplabel(['C test, alpha_C=' num2str(10^alpha_Cs(i))], 't');
    
    R = inten_test(1:32,1:end-2) - langmuir_func(A, B, C_test);
    test_errors_fro(i) = sum(sum(R.^2)) /  test_size;
    %test_errors_l1(i)  = sum(sum(abs(R))) / test_size;
    %test_errors_weighted(i) = sum(sum(loss_asymmetric(R) .* inten_test_sliced{2} )) / test_size;
    
    if ((err_best-test_errors_fro(i))/test_errors_fro(i)>0.001)||...
            (err_best-test_errors_fro(i)>0 && reg_best>=10^alpha_Cs(i))
        Abest = A;
        Bbest = B;
        Cbest = C;
        Abestsliced = A_sliced;
        Bbestsliced = B_sliced;
        Ctestbest = C_test;
        reg_best = 10^alpha_Cs(i);
        err_best = test_errors_fro(i);
    end
    
    if i>1
        if (test_errors_fro(i)-test_errors_fro(i-1))/test_errors_fro(i-1)>0.01
            goup = ~goup;
        end
        if goup
            prevUps = alpha_Cs(alpha_Cs>alpha_Cs(i));
            if isempty(prevUps)
                alpha_Cs(i+1) = alpha_Cs(i) + 5;
            else
                if (min(prevUps)==alpha_Cs(i)+1)
                    ind = find(alpha_Cs == min(prevUps));
                    if (test_errors_fro(ind)-test_errors_fro(i))/test_errors_fro(ind) >0.01%(alpha_Cs(ind-1)<alpha_Cs(ind))
                        goon = false;
                    else
                        alpha_Cs(i+1) = min(setdiff(alpha_Cs(i)+1:max(alpha_Cs(2:end)), prevUps));
                    end
                else
                    alpha_Cs(i+1) = ceil((alpha_Cs(i)+min(prevUps))/2);
                end
            end
        else
            if abs(test_errors_fro(i)-test_errors_fro(1))/test_errors_fro(1)<0.001
                goon = false;
            else
                prevDowns = alpha_Cs(alpha_Cs<alpha_Cs(i));
                prevDowns = prevDowns(~isinf(prevDowns));
                if isempty(prevDowns)                    
                    alpha_Cs(i+1) = alpha_Cs(i) - 5;
                else
                    if (max(prevDowns)==alpha_Cs(i)-1)
                        ind = find(alpha_Cs == max(prevDowns));
                        if (test_errors_fro(ind)-test_errors_fro(i))/test_errors_fro(ind)>0.01%(alpha_Cs(ind-1)>alpha_Cs(ind))
                            goon = false;
                        else
                            alpha_Cs(i+1) = max(setdiff(min(alpha_Cs(2:end)):alpha_Cs(i)-1, prevDowns));
                        end
                    else
                        alpha_Cs(i+1) = floor((alpha_Cs(i)+max(prevDowns))/2);
                    end
                end
            end
        end
    end
    
    i=i+1;
end
% 
%  for i=1:length(alpha_Cs)
%      alpha_Cs(i)
%      [A, B, C] =  nonlinear_alpha_beta_LO_reg(inten_train_sliced{1}, alpha, beta, maxIterCnt, eps, alpha_Cs(i), 1);
%      C_test = nonlinear_alpha_beta_fixedAB(inten_test_sliced{1}, A, B, alpha, beta, maxIterCnt, eps, alpha_Cs(i), 1);
% 
%      R = inten_test_sliced{1} - langmuir_func(A, B, C_test);
%      test_errors_fro(i) = sum(sum(R.^2)) /  size(inten_test_sliced{1},2);
% %     test_errors_l1(i)  = sum(sum(abs(R))) / size(inten_test_sliced{2},2);
% %     test_errors_weighted(i) = sum(sum(loss_asymmetric(R) .* inten_test_sliced{2} )) / size(inten_test_sliced{2},2);
%  end

