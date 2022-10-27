function dex = ten_exp_d1(v,v_d)
    % ten_exp_d gives the (d(exp(theta))/dS)
    
    e = norm(v);
    
    if e == 0
       dex1 = 1*(skew_ten(v_d)) + ((1/2)*(skew_ten(v_d)*skew_ten(v) + skew_ten(v)*skew_ten(v_d)));
    else
       dex1 = ((sin(e)/e)*(skew_ten(v_d))) + (((1-cos(e))/(e^2))*(skew_ten(v_d)*skew_ten(v) + skew_ten(v)*skew_ten(v_d))) + (((e*cos(e) - sin(e))/(e^2))*(dot(v,v_d)/e)*skew_ten(v)) + (((e*sin(e) - 2 + 2*cos(e))/(e^3))*(dot(v,v_d)/e)*(skew_ten(v)^2));                                                          
    end
    
    dex = dex1;%*dex1';%ten_exp(-v);

end