function dex = ten_exp_d(v,v_d)
    % ten_exp_d gives the (d(exp(theta))/dS)*exp(-theta)
    h = 0.001;
    if norm(v) == 0
        e = v/h;
    else
        e = v/norm(v);
    end
    v_bar = tan(norm(v)/2)*e;
    theta_bar = skew_ten(v_bar);
    
    if norm(v) == 0
        v_bar_d = ((1/2)*1)*(v_d - ((1-1)*(e'*v_d)*e));
    else
        v_bar_d = ((1/2)*(tan(norm(v)/2))/((1/2)*norm(v)))*(v_d - ((1-(norm(v)/sin(norm(v))))*(e'*v_d)*e));
    end
    theta_bar_d = skew_ten(v_bar_d);
    dex = (2/(1+norm(v_bar)^2));%*(theta_bar_d + theta_bar_d*theta_bar - theta_bar*theta_bar_d);

end