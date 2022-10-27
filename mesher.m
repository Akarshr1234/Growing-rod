function z_cor = mesher(r,n,l)
    a = l*(1-r)/(1-r^n);
    z_cor = zeros(1,n+1);
    for node = 2:n+1
        z_cor(node) = a*(1-r^(node-1))/(1-r);
    end
end