function min_point = ternary_search(func, left, right, eps)
    phi = 1.61803399;
    
    while (right - left > eps)
        %m1 = left + (right - left) / 3;
        %m2 = left + 2 * (right - left) / 3;
        m1 = right - (right - left) / phi;
        m2 = left + (right - left) / phi;
        
        if (func(m1) < func(m2))
            right = m2;
        else
            left = m1;
        end
    end
    
    min_point = (left + right) / 2;
end