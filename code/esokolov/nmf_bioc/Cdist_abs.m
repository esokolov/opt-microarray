function d = Cdist_abs(C1,C2)
    good_ind = ~isnan(C1) & ~isnan(C2) & ~isinf(C1) & ~isinf(C2);
    C1 = C1(good_ind);
    C2 = C2(good_ind);
    if size(C1,1)==1 
        C1=C1'; 
    end
    if size(C2,1)==1 
        C2=C2'; 
    end
    C = [C1 C2];
    cov_matr = C' * C;
    if (sum(sum(isinf(cov_matr))) == 0)
        [V,~] = eig(C' * C);
        d = mean(abs(C * V(:,1)));
    else
        d = 1e100;
    end
end