function vari = NR_iter(res_fun,vari)

    global nnodes;
    global l;
    global Force;
    global len_ini;
    global pos_ini;
    global xhi_ini;
    global body_Force;
    niter = 15;                     % number of iterations
    % Newton Raphson Iteration
%     pos = vari(1:3*nnodes);
%     xhi = vari(3*nnodes+1:end);
%     plot(pos(3:3:end),pos(1:3:end),'LineWidth',1,'color','b','LineStyle','--');
%     pause(1);
%     len = len_finder(pos,nel);
     for iter = 1:niter
               %% getting solution of variation
                [G,Tangent] = feval(res_fun,vari);
                Delta_var = - Tangent\(G-(Force+body_Force));
%                 Delta_pos = Delta_var(1:end/2);
%                 Delta_xhi = Delta_var((end/2)+1:end);
                %% update
                vari = vari + Delta_var;
                pos = vari(1:3*nnodes);
                xhi = vari(3*nnodes+1:end);
                if norm(G-(Force+body_Force)) <= 10^-10
                    break;
                end

     end 
%                 plot(pos(3:3:end),pos(1:3:end),'LineWidth',1,'color','b','LineStyle','--');
                norm(G-(Force+body_Force))
                iter
%                 pause(1);
%                 clf;

    %  norm(G-Force)
    end
    