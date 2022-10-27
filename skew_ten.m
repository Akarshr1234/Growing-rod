function sk = skew_ten(vec)
% skew gives the skew matrix corresponding to size 3 column vector
    if length(vec) ~= 3
        error('skew function only works with a vector of size 3');
    end
    sk = zeros(3,3);
    sk(1,2) = -vec(3);
    sk(1,3) = vec(2);
    sk(2,3) = -vec(1);
    
    sk = sk - sk';
end