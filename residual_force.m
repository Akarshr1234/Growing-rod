function [G,Tangent,kappa_all_gp] = residual_force(variables)

    global l;            % total length
    global nel;          % number of locals
    global npe;          % number of nodes per local
    global nnodes;       % total number of nodes
    global ndof;         % total number of positional degree of freedom
    global ngp;          % number of gauss points to be used

    global n_bar;        % initialisation of body force
    global m_bar;        % initialisation of body moment

    global bcs;          % bounadary condition dof
    global bc_val;       % boundary condition values
    global C;            % material parameters
    global Be_stiffness; % EI of each elements
    global pos_ini;      % undeformed configuration position
    global xhi_ini;      % undeformed configuration rotation

    global Force;        % natural boundary condition
    global body_Force;   % body force vector
    global L;            % gather operator

    pos = variables(1:end/2);
    xhi = variables((end/2)+1:end);

    len = len_finder(pos_ini,nel);                  % length of each local
    [gps,w] = gauss_legendre(ngp);              % gauss points and weights

    G = zeros(nnodes*6,1);                      % global residual vector initialisation
    body_Force = zeros(nnodes*6,1);             % body forces
    Tangent = zeros(nnodes*6,nnodes*6);         % global tangent matrix initialisation
    kappa_all_gp = zeros(nel,ngp);              % kappa at gauss points for growth
    Pi_axial = zeros(3*npe,1);                  % local nodal xhi initialisation
    Pi_axial_ini = zeros(3*npe,1);
    for el = 1:nel 
       pos_el = L(:,:,el)*pos;                  % local nodal position
       pos_el_ini = L(:,:,el)*pos_ini;
       xie = zeros(npe*3,npe*3);                % local nodal rotation
       C_el = C;
       C_el(5,5) = Be_stiffness(el);
       for node = 1:npe
           p = 3*(((npe-1)*(el-1)+node)-1);                             % index
           Pi_axial(3*(node-1)+1 : 3*(node-1)+3) = xhi(p+1:p+3);       % local nodal xhi definition
           Pi_axial_ini(3*(node-1)+1 : 3*(node-1)+3) = xhi_ini(p+1:p+3);
       end

       g = zeros(npe*3+npe*3,1);                                       % local residual vector
       body_force = zeros(npe*3+npe*3,1);
       tangent_ele = zeros(npe*3+npe*3,npe*3+npe*3);                              % local element tangent stiffness
       tangent_geo = zeros(npe*3+npe*3,npe*3+npe*3);                              % local geometric tangent stiffness
       n_bar_el = L(:,:,el)*n_bar;                               % body force of element
       m_bar_el = L(:,:,el)*m_bar;                               % body moment of element


       % gauss integration
       for gp = 1:ngp
          N_shape = shape_fun(gps(gp),npe-1);                           % shape function at gauss point
          N_shape_der = shape_fun_der(gps(gp),npe-1)*(2/len(el));       % shape function derivative at gauss point
          shape_der = zeros(3,npe*3);                                   % summation of N_shape_der*eye(3) matrix initialisation
          shape_ele = zeros(3,npe*3);                                   % summation of N_shape*eye(3) matrix initialisation
          for node = 1:npe
                shape_der(:,3*(node-1)+1:3*node) = diag(N_shape_der(node)*[1 ,1 ,1]);       % summation of N_shape_der*eye(3) matrix definition
                shape_ele(:,3*(node-1)+1:3*node) = diag(N_shape(node)*[1, 1, 1]);           % summation of N_shape*eye(3) matrix definition
          end       
          der_phi = shape_der*pos_el;                          % derivative of position at gauss point
          der_phi_ini = shape_der*pos_el_ini;
          der_xi = shape_der*Pi_axial;                          % derivative of rotation at gauss point
          Lam_ax_gp = shape_ele*Pi_axial;                       % xhi at gauss point
          Lam_gp = ten_exp(Lam_ax_gp);                          % lamda at gauss point

          der_xi_ini = shape_der*Pi_axial_ini;
          Lam_ax_gp_ini = shape_ele*Pi_axial_ini; 
          Lam_gp_ini = ten_exp(Lam_ax_gp_ini);

          Omega_gp_ini =  ten_exp_d1(Lam_ax_gp_ini,der_xi_ini)*ten_exp(Lam_ax_gp_ini)' ;
          Omega_gp = ten_exp_d1(Lam_ax_gp,der_xi)*ten_exp(Lam_ax_gp)'  +  Lam_gp*Omega_gp_ini*Lam_gp';      % omega at gauss point
          kappa_gp = Lam_gp_ini'*Lam_gp'*(axial(Omega_gp)); %                             % kappa at gauss point
          kappa_all_gp(el,gp) = kappa_gp(2,1);
          Gamma_gp = Lam_gp_ini'*(Lam_gp'*der_phi) - [0;0;1]; %                           % gamma at gauss point
          Reac = 2*(1/2)*C_el*[Gamma_gp - (Lam_gp_ini'*der_phi_ini - [0;0;1]) ;kappa_gp - (Lam_gp_ini'*(axial(Omega_gp_ini))) ];                             % reactions at gauss point 
          N_gp = Reac(1:3);                                                 % nodal force at gauss piont
          M_gp = Reac(4:6);                                                 % nodal moments at gauss point

          phid_cr = zeros(3,npe*3);             % summation of N_shape * cross with derivative of position at gauss point initialisation
          for node = 1:npe
             phid_cr(:,3*(node-1)+1:3*node) = N_shape(node)*skew_ten(der_phi);      % summation of N_shape * cross with derivative of position at gauss point definition
          end

          % shape function matrix for position and rotation
          shape_fun_ele = [shape_ele, zeros(size(shape_ele));                   
                            zeros(size(shape_ele)),shape_ele];



          shape_R = zeros(3,npe*3);

          for node = 1:npe
              shape_R(:,3*(node-1)+1:3*(node-1)+3) = N_shape(node)*eye(3);
          end

          n_gp = Lam_gp*Lam_gp_ini*N_gp;                           % body force at gauss point
          m_gp = Lam_gp*Lam_gp_ini*M_gp;                           % body moment at gauss point
          Nn_cr = zeros(3,npe*3);                       % N_shape[n cross] operator at gauss point initialisation
          Nm_cr = zeros(3,npe*3);                       % N_shape[m cross] operator at gauss point initialisation
          Ndn_cr = zeros(3,npe*3);                      % N_shape_der[n cross] operator at gauss point initialisation
          Nnp_cr = zeros(3,npe*3);                      % N_shape([n dydic pos_der] - (n.pos_der)*eye(3)) operator at gauss point initialisation
          for node = 1:npe
             Nn_cr(:,3*(node-1)+1:3*node) = N_shape(node)*skew_ten(n_gp);       % N_shape[n cross] operator at gauss point definition
             Nm_cr(:,3*(node-1)+1:3*node) = N_shape(node)*skew_ten(m_gp);       % N_shape[m cross] operator at gauss point definition
             Ndn_cr(:,3*(node-1)+1:3*node) = N_shape_der(node)*skew_ten(n_gp);  % N_shape_der[n cross] operator at gauss point definition
             Nnp_cr(:,3*(node-1)+1:3*node) = N_shape(node)*(dydic(n_gp,der_phi)-(n_gp'*(der_phi)*eye(3))); % N_shape([n dydic pos_der] - (n.pos_der)*eye(3)) operator at gauss point definition
          end

          % xi at gauss point
          xi_gp = [shape_der,phid_cr;
                   zeros(size(shape_der)),shape_der];
          % psi at gauss point
          psi_gp = [shape_der, zeros(size(shape_der));
                    zeros(size(shape_der)), shape_der;
                    zeros(size(shape_der)), shape_der];
          % be at gauss point
          Be_gp = [zeros(size(Nn_cr)), -Nn_cr;
                   zeros(size(Nm_cr)), -Nm_cr;
                   Ndn_cr, Nnp_cr];
          % pi at gauss point
          Pi_gp = [Lam_gp*Lam_gp_ini, zeros(size(Lam_gp));
                   zeros(size(Lam_gp)), Lam_gp*Lam_gp_ini];
          g = g + (len(el)/2)*w(gp)*((xi_gp'*[n_gp;m_gp]));          % residual at gauss point
          body_force = body_force + (len(el)/2)*w(gp)*((shape_fun_ele'*shape_fun_ele)*[n_bar_el;m_bar_el]); 
          tangent_ele = tangent_ele + (len(el)/2)*w(gp)*(xi_gp'*Pi_gp*C_el*Pi_gp'*xi_gp);                                          % element tangent stiffness at gauss point
          tangent_geo = tangent_geo + (len(el)/2)*w(gp)*(psi_gp'*Be_gp);                                                        % geometric tangent stiffness at gauss point
       end

       tangent_stiff = tangent_ele + tangent_geo;                       % total tangent stiffness at gauss point

       %% Assembly
       % gather operator for position and rotation
       L_whole = [L(:,:,el), zeros(size(L(:,:,el)));
                  zeros(size(L(:,:,el))), L(:,:,el)];
       % global tangent stiffness matrix 
       Tangent = Tangent + L_whole'*tangent_stiff*L_whole;
       % global residual vector
       G = G + L_whole'*g;
       body_Force = body_Force + L_whole'*body_force;
    end

    %% imposing essential boundary condition
    for bc = 1:length(bcs)
        G(bcs(bc)) = 0;
        body_Force(bcs(bc)) = 0;
        Tangent(bcs(bc),:) = zeros(size(Tangent(bcs(bc),:)));
        Tangent(:,bcs(bc)) = zeros(size(Tangent(:,bcs(bc))));
        Tangent(bcs(bc),bcs(bc)) = 1;
    end
end



        
