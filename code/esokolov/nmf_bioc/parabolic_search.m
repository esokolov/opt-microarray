function min_point = parabolic_search(func, left, right, eps)
    m = (left + right) / 2;
    
    while (right - left > eps)
        if (left == m)
            m = left + eps;
        end
        if (right == m)
            m = right - eps;
        end
        
        f_left = func(left);
        f_m = func(m);
        f_right = func(right);
        
        A = [left^2 left 1;
             m^2 m 1;
             right^2 right 1 ];
        b = [func(left) func(m) func(right)]';
        coefs = A \ b;
        
        new_point = max(-coefs(2) / (2 * coefs(1)), 0);
        f_new = func(new_point);
        
        if (new_point < m)
            if (f_new < f_m)
                right = m;
                m = new_point;                
            else
                left = new_point;
            end
        else
            if (f_new < f_m)
                left = m;
                m = new_point;
            else
                right = new_point;
            end
        end
    end
    
    min_point = (left + right) / 2;
end