function var = arc_length(res_fun,var)
    psi = 1;
    Delta_l = 1;
    tol = 10^-4;
    niter = 20;
    lambda_max = 1;
    lambda = 0.05;
    global nnodes;
    global Force;
    global len_ini;
    global pos_ini;
    global xhi_ini;
    global l;
    global body_Force;
    
    while lambda<lambda_max
        [G,Tangent] = feval(res_fun,var);
        Delta_var = zeros(6*nnodes,1);
        Delta_lambda = 0;
        delta_var_bar = zeros(6*nnodes,1);
        delta_var_t = Tangent\(Force + body_Force);

        alpha_1 = delta_var_t'*delta_var_t + psi^2*((Force + body_Force)'*(Force + body_Force));
        alpha_2 = 2*(Delta_var + delta_var_bar)'*delta_var_t + 2*psi^2*Delta_lambda*((Force + body_Force)'*(Force + body_Force));
        alpha_3 = (Delta_var + delta_var_bar)'*(Delta_var + delta_var_bar) + (psi*Delta_lambda)^2*((Force + body_Force)'*(Force + body_Force)) - Delta_l^2;
        if alpha_1 == 0
            delta_lambda = -alpha_3/alpha_2;
        else
            if det(Tangent)>0
                delta_lambda = real(-alpha_2 + sqrt(alpha_2^2 - 4*alpha_3*alpha_1))/(2*alpha_1);
            else
                delta_lambda = real(-alpha_2 - sqrt(alpha_2^2 - 4*alpha_3*alpha_1))/(2*alpha_1);
            end
        end
        delta_var = delta_var_bar + delta_lambda*delta_var_t;

        var = var + delta_var;
        lambda = lambda+delta_lambda;

        [G,Tangent] = feval(res_fun,var);
        res = G - lambda*(Force + body_Force);
        if norm(res) > tol
            Delta_var = delta_var;
            Delta_lambda = delta_lambda;
            for iter =1:niter
                delta_var_bar = -Tangent\(G-lambda*(Force + body_Force));
                delta_var_t = Tangent\(Force + body_Force);
                p_Delta_labmda = Delta_lambda;
                p_Delta_var = Delta_var;
                

                alpha_1 = delta_var_t'*delta_var_t + psi^2*((Force + body_Force)'*(Force + body_Force));
                alpha_2 = 2*(Delta_var + delta_var_bar)'*delta_var_t + 2*psi^2*Delta_lambda*((Force + body_Force)'*(Force + body_Force));
                alpha_3 = (Delta_var + delta_var_bar)'*(Delta_var + delta_var_bar) + (psi*Delta_lambda)^2*((Force + body_Force)'*(Force + body_Force)) - Delta_l^2;
                

                delta_lambda_1 = real(-alpha_2 + sqrt(alpha_2^2 - 4*alpha_3*alpha_1))/(2*alpha_1);
                delta_lambda_2 = real(-alpha_2 - sqrt(alpha_2^2 - 4*alpha_3*alpha_1))/(2*alpha_1);
                delta_var_1 = delta_var_bar + delta_lambda_1*delta_var_t;
                delta_var_2 = delta_var_bar + delta_lambda_2*delta_var_t;
                
                arc1 = (Delta_var + delta_var_1)'*Delta_var + psi^2*(Delta_lambda + delta_lambda_1)*Delta_lambda*(Force + body_Force)'*(Force + body_Force);
                arc2 = (Delta_var + delta_var_2)'*Delta_var + psi^2*(Delta_lambda + delta_lambda_2)*Delta_lambda*(Force + body_Force)'*(Force + body_Force);
                
                if arc1 > arc2
                    delta_var = delta_var_1;
                    delta_lambda = delta_lambda_1;
                else
                    delta_var = delta_var_2;
                    delta_lambda = delta_lambda_2;
                end
                
                
                Delta_var = delta_var;
                Delta_lambda = delta_lambda;
                var = var + Delta_var;
                lambda = lambda + Delta_lambda;
                if lambda > lambda_max
                    lambda = lambda_max;
                    return;
                end 
                [G,Tangent] = feval(res_fun,var);
                res = norm(norm(G - lambda*(Force + body_Force)));
                if res<tol
                    break;
                end

            end
        end
        norm(res)
        lambda
        iter
        pos = var(1:end/2);
        xhi = var((end/2)+1:end);
%         pos(2:3:end)
%         pause(0.2);
        plot(pos_ini(3:3:end),pos_ini(1:3:end),'r');
        hold on;
        plotter(xhi_ini,pos_ini,l/10,'r');
        plot(pos(3:3:end),pos(1:3:end),'b');
        xlim([-sum(len_ini),sum(len_ini)]);
        ylim([-sum(len_ini),sum(len_ini)]);
        plotter(xhi+xhi_ini,pos,l/10,'b');
        drawnow();
%         pause(0.01);
        clf;
               
%         axis equal;
%         axis equal;
    end
end 


