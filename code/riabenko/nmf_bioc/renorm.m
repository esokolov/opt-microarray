function [A_new C_new isbreak] = renorm(I, A, C, maxIterCnt, eps)
%на вход подаются матрицы интенсивностей, векторы аффинитивностей и
%концентраций, причём sum(abs(A))==1; maxIterCnt - число проходов по
%пробсету, eps - максимальное изменение концентраций от перенормировки, при
%котором процесс всё ещё продолжается, isbreak - индикатор сходимости
%концентраций.
A_new = A;
C_new = C;
isbreak = false;
for currIter = 1:maxIterCnt
    for i=1:length(A)
        A_loo = A_new(setdiff(1:end,i)) / sum(abs(A_new(setdiff(1:end,i))));
        C_loo = nmf_smart_fixedA(I(setdiff(1:end,i),:), A_loo,100,10^-6);
        A_new(i) = nmf_smart_fixedC(I(i,:),C_loo,100,10^-6);
        A_new = A_new / sum(abs(A_new));
    end
    C_prev = C_new;
    C_new = nmf_smart_fixedA(I, A_new,100,10^-6);
    if sum((C_prev-C_new).^2)<eps
        isbreak = true;
        break;
    end
end

