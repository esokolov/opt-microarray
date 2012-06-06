function point = generateFromCdf(x, cdf)
    r = rand(1);
    b = find(r < cdf, 1, 'first');
    a = b - 1;
    if (a == 0)
        ya = 0;
        xa = 0;
    else
        ya = cdf(a);
        xa = x(a);
    end
    xb = x(b);
    yb = cdf(b);
    point = xa + ((xb - xa) / (yb - ya)) * (r - ya);
end