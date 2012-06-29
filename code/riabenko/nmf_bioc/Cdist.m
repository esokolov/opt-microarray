function d = Cdist(C1,C2)
if size(C1,1)==1 
    C1=C1'; 
end
if size(C2,1)==1 
    C2=C2'; 
end
C = [C1 C2];
[V,~] = eig(C' * C);
d = sum((C * V(:,1)).^2);