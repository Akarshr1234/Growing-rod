function [pos_ini,xhi_ini] = update_curvature(xhi_ini,pos_ini,xhi,pos)
    global nel;
    global ngp;
    global npe;
    global nnodes;
    global Lo_dot;
    global beta_growth; 
    global gamma_growth;
    global lg_growth;
    global dt;
    global Be_stiffness;
    global delta_B_stiffness;
    global tau;
    global age_vector;

    [~,~,kappa_gp_all] = residual_force([pos;xhi]);
    len_ini = len_finder(pos_ini,nel);                  % length of each local
    [gps,w] = gauss_legendre(ngp);              % gauss points and weights
    kappa_ini_gp = zeros(nel,ngp);
    Pi_axial_ini = zeros(3*npe,1);
    Pi_axial = zeros(3*npe,1);
    xhi_gp_all = zeros(nel,ngp); 
    for el = 1:nel
        for node = 1:npe
           p = 3*(((npe-1)*(el-1)+node)-1);                             % index
           Pi_axial_ini(3*(node-1)+1 : 3*(node-1)+3) = xhi_ini(p+1:p+3);
           Pi_axial(3*(node-1)+1 : 3*(node-1)+3) = xhi(p+1:p+3);
       end

        for gp = 1:ngp
            N_shape = shape_fun(gps(gp),npe-1);                           % shape function at gauss point
            N_shape_der = shape_fun_der(gps(gp),npe-1)*(2/len_ini(el));       % shape function derivative at gauss point
            shape_der = zeros(3,npe*3);                                   % summation of N_shape_der*eye(3) matrix initialisation
            shape_ele = zeros(3,npe*3);                                   % summation of N_shape*eye(3) matrix initialisation
            for node = 1:npe
               shape_der(:,3*(node-1)+1:3*node) = diag(N_shape_der(node)*[1 ,1 ,1]);       % summation of N_shape_der*eye(3) matrix definition
               shape_ele(:,3*(node-1)+1:3*node) = diag(N_shape(node)*[1, 1, 1]);           % summation of N_shape*eye(3) matrix definition
            end
            der_xi_ini = shape_der*Pi_axial_ini;
            Lam_ax_gp_ini = shape_ele*Pi_axial_ini; 
            Lam_gp_ini = ten_exp(Lam_ax_gp_ini);
            Omega_gp_ini =  ten_exp_d1(Lam_ax_gp_ini,der_xi_ini)*ten_exp(Lam_ax_gp_ini)';
            xhi_gp = shape_ele*(Pi_axial+Pi_axial_ini);
            kappa_gp = Lam_gp_ini'*(axial(Omega_gp_ini));
            kappa_ini_gp(el,gp) = kappa_gp(2);
            xhi_gp_all(el,gp) = xhi_gp(2);
        end
    end
%     kappa_ini_gp
    kappa_ini_gp = kappa_ini_gp + (Lo_dot*dt*((-beta_growth*sin(xhi_gp_all))+(-gamma_growth*(kappa_gp_all))));
%     kappa_ini_gp
%     pause(5);
    int_kappa_ini = zeros(nel,1);
    gr_in = growth_index(len_ini,lg_growth);
    for el = 1:nel
        for gp = 1:ngp
            int_kappa_ini(el) = int_kappa_ini(el) + w(gp)*(len_ini(el)/2)*kappa_ini_gp(el,gp);
        end
    end
    int_kappa_ini_req = int_kappa_ini(gr_in:end);
    new_xhi_ini = xhi_ini(3*(gr_in-1)+2)+[0;cumsum(int_kappa_ini_req)];
    for node = gr_in:nnodes
        xhi_ini(3*(node-1)+2) = new_xhi_ini(node - gr_in +1);
    end
    for el = gr_in:nel
       len_ini(el) = len_ini(el)*(1 + Lo_dot*dt);  
    end
%         
    for node = (gr_in+1):nnodes
        pos_ini(3*(node-1)+1:3*(node-1)+3) = pos_ini(3*(node-2)+1:3*(node-2)+3) + (ten_exp(xhi_ini(3*(node-2)+1:3*(node-2)+3)))*[0;0;len_ini(node-1)];
    end
    if gr_in>1
        age_vector(1:gr_in-1) = age_vector(1:gr_in-1) + 1;
    end
    for el = 1:gr_in-1
        Be_stiffness(el) = Be_stiffness(el) + delta_B_stiffness*(exp((-age_vector(el)*dt)/tau) - exp(-((age_vector(el)+1)*dt)/tau));
    end

end