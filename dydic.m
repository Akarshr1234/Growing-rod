function ten = dydic(a,b)
    if length(a) ~= length(b)
        error('Size of both inputs are not same');
    end
    
    ten = zeros(length(a),length(b));
    for k = 1:length(a)
        for p = 1:length(b)
            ten(k,p) = a(k)*b(p);
        end
    end
    
end

