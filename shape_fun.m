function Nx = shape_fun (x,xo)
    x0 = -1:(1-(-1))/xo:1;

    
    Nx = zeros(1, length(x0));

    
    for i = 1:length(x0)
        denom = 1;
        numer = 1;
        for j = 1:length(x0)
            if j == i 
                continue;
            end
            denom = denom*(x0(j) - x0(i));
            numer = numer*(x0(j) - x);
        end
        Nx(i) = numer/denom;
    end 
end




