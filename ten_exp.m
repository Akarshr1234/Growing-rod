function ex = ten_exp(v)
    % find exponential of a skew tensor corresponding to a 3X1 vector.
    h = 0.0001;
    if norm(v) == 0
        e = v/h;
    else
        e = v/norm(v);
    end
    v_bar = tan(norm(v)/2)*e;
    theta_bar = skew_ten(v_bar);
    
    ex = eye(3) + (2/((1+(norm(v_bar))^2)))*(theta_bar + theta_bar^2);

end