function Nx = shape_fun_der (x,xo)
    x0 = -1:(1-(-1))/xo:1;

    
    Nx = zeros(1, length(x0));
    
    for i = 1:length(x0)
       denom = 1;
       for j = 1:length(x0)
           if j == i
               continue;
           end
           denom = denom*(x0(j) - x0(i));
       end
       
       term = 0;
       
       for j = 1:length(x0)
           if j == i
               continue;
           end
          numer = -1;
          for k = 1:length(x0)
              if k == i || k == j
                  continue;
              end
              numer = numer*(x0(k) - x);
          end
          term = term + (numer/denom);
       end
       Nx(i) = term;
    end
    
    
end




