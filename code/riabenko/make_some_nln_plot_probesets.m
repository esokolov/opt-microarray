as = [6,  6, 8, 9];
bs = [11, 4, 5, 4];
for i=1:length(as)
    C = nonlinear_alpha_beta_LO_reg_fixedAB(inten_full_sliced{2}+1, A_sliced_all{as(i),bs(i)}{2}, B_sliced_all{as(i),bs(i)}{2}, ...
        alpha_range(as(i)), beta_range(bs(i)), maxIterCnt, eps, reg_best(as(i),bs(i)), 1);
    nln_plot_probeset(inten_full_sliced{2}+1,A_sliced_all{as(i),bs(i)}{2}, B_sliced_all{as(i),bs(i)}{2},C,0)
end