function d = Cdist(C1,C2)
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
        d = sum((C * V(:,1)).^2);
    else
        d = 1e100;
    end
end