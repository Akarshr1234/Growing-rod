l = 10;                         % total length
nel = 5;                       % number of locals
npe = 2;                        % number of nodes per local
nnodes = (npe-1)*nel + 1;       % total number of nodes
ndof = 3* nnodes;               % total number of positional degree of freedom
ngp = 1;                        % number of gauss points to be used
nor = [];                       % norm of residual tracking
G_iter = [];                    % residual tracking
vari_iter = [];                 % variation update tracking
pos_iter = [];                  % position update tracking
% n_bar and m_bar value
n_bar_c = [0; 0; 0];            % constant body force vector
m_bar_c = [0; 0; 0];            % constant body moment vector
n_bar = zeros(nnodes*3,1);      % initialisation of body force
m_bar = zeros(nnodes*3,1);      % initialisation of body moment

%% Essential Boundary Conditions
bcs = [1,2,3,nnodes*3+1,nnodes*3+2,nnodes*3+3];%,6*nnodes,6*nnodes-1,6*nnodes-2,nnodes*3-2,nnodes*3-1]; % bounadary condition dof
bc_val = [0,0,0,0,0,0,0,0,0,0,0,0]; % boundary condition values

% body forces and body moment definition
for node = 1:nnodes
    n_bar(((node-1)*3+1):((node-1)*3+3)) = n_bar_c;
    m_bar(((node-1)*3+1):((node-1)*3+3)) = m_bar_c;
end


%% Problem Setup
L = zeros(npe*3,ndof,nel);                      % gather operator ititialisation
% C = diag([10^7 10^8 10^6 10^4 10^5 10^3]);
% C = diag((1/pi^2)*ones(6,1));                        % material parameters
C = diag(100*ones(6,1));     
% C(1,1) = 100*10^5;
% C(2,2) = C(1,1);
% C(3,3) = 100*10^5;

% C = diag([1.61538 *10^8,1.61538 *10^8,3.5 *10^7,3.5 *10^7,3.5 *10^7,3.5 *10^7]);
% C = diag([3600*10^10,3600*10^10,7200*10^10,7200*0.5,7200*0.5,7200]);


% gather operator definition
for el = 1:nel
    L(:,(el-1)*(3*(npe-1)) + 1 :((el-1)*(3*(npe-1))) + npe*3,el) = eye(npe*3);
end

% initial position discription
posz = linspace(0,l,nnodes)';
% posz = mesh(0.9,nel,l);
pos = zeros(nnodes*3,1);
for node = 1:nnodes
   pos(3*(node)) = posz(node); 
end

% initial rotation description

xhi = zeros(nnodes*3,1);
xhi_ini = zeros(nnodes*3,1);

% initial variation guess
uo = zeros(nnodes*3,1);
nu = zeros(nnodes*3,1);

len = len_finder(pos,nel);                  % length of each local
len_ini = len;
[gps,w] = gauss_legendre(ngp);              % gauss points and weights
    
dt = 1;
time = 1;
tsteps = time/dt;
timeins = linspace(0,time,tsteps);
xlim([-sum(len_ini),sum(len_ini)]);
ylim([-sum(len_ini),sum(len_ini)]);
drawnow;
pause(1);


len_iter = sum(len);
for t = 1:length(timeins)
    
    clc;
    fprintf('Time = %f \n',timeins(t));
    fprintf('Timestep = %d \n',t);
    xhi_ini = xhi + xhi_ini;
    xhi = zeros(nnodes*3,1);
    pos_ini = pos;
    len = len_finder(pos,nel);                  % length of each local
    niter = 100;                     % number of iterations
    % Newton Raphson Iteration
    for iter = 1:niter

    %     iter

        G = zeros(nnodes*6,1);                      % global residual vector initialisation  
        Tangent = zeros(nnodes*6,nnodes*6);         % global tangent matrix initialisation
        Pi_axial = zeros(3*npe,1);                  % local nodal xhi initialisation
        Pi_axial_ini = zeros(3*npe,1);
        for el = 1:nel 
           pos_el = L(:,:,el)*pos;                  % local nodal position
           pos_el_ini = L(:,:,el)*pos_ini;
           xie = zeros(npe*3,npe*3);                % local nodal rotation
           for node = 1:npe
               p = 3*(((npe-1)*(el-1)+node)-1);                             % index
               Pi_axial(3*(node-1)+1 : 3*(node-1)+3) = xhi(p+1:p+3);       % local nodal xhi definition
               Pi_axial_ini(3*(node-1)+1 : 3*(node-1)+3) = xhi_ini(p+1:p+3);
           end

           g = zeros(npe*3+npe*3,1);                                       % local residual vector
           tangent_ele = zeros(npe*3+npe*3,npe*3+npe*3);                              % local element tangent stiffness
           tangent_geo = zeros(npe*3+npe*3,npe*3+npe*3);                              % local geometric tangent stiffness
           n_bar_el = L(:,:,el)*n_bar;                               % body force of element
           m_bar_el = L(:,:,el)*m_bar;                               % body moment of element


           % gauss intigration
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
              Gamma_gp = Lam_gp_ini'*(Lam_gp'*der_phi) - [0;0;1]; %                           % gamma at gauss point
              Reac = 2*(1/2)*C*[Gamma_gp - (Lam_gp_ini'*der_phi_ini - [0;0;1]) ;kappa_gp - (Lam_gp_ini'*(axial(Omega_gp_ini))) ];                             % reactions at gauss point 
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
              g = g + (len(el)/2)*w(gp)*((xi_gp'*[n_gp;m_gp]) - (shape_fun_ele'*shape_fun_ele)*[n_bar_el;m_bar_el]);          % residual at gauss point
              tangent_ele = tangent_ele + (len(el)/2)*w(gp)*(xi_gp'*Pi_gp*C*Pi_gp'*xi_gp);                                          % element tangent stiffness at gauss point
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

        end

        %% imposing essential boundary condition
        for bc = 1:length(bcs)
            G(bcs(bc)) = 0;
            Tangent(bcs(bc),:) = zeros(size(Tangent(bcs(bc),:)));
            Tangent(:,bcs(bc)) = zeros(size(Tangent(:,bcs(bc))));
            Tangent(bcs(bc),bcs(bc)) = 1;
        end
        %% imposing natural boundary condition
    %     G(3*(((nnodes + 1)/2))-2) = G(3*(((nnodes + 1)/2))-2) - 0.00000002;
%         G(3*(((nnodes + 1)/2))) = G(3*(((nnodes + 1)/2))) - 2;
%         G(3*(nnodes)-2) = G(3*(nnodes)-2) - 0.02;
    % %       temp = G(3*nnodes)
%           G(3*nnodes) = G(3*nnodes) + 1*pi^2*C(4,4)/l;
%           G(3*nnodes) = G(3*nnodes) - 0.641713*10^-1;
%           G(3*nnodes) = G(3*nnodes) - 0.002;

%           G(3*nnodes) = G(3*nnodes) + 1;

    %       G(3*((nnodes+1)/2)) = G(3*((nnodes+1)/2)) - 1;

    %         G(6*nnodes) = G(6*nnodes) - 100;

%         G(6*nnodes-1) = G(6*nnodes-1)-((C(4,4)*2*pi/l)/(2));
        %% getting solution of variation
        variables = -Tangent\G;
        uo = variables(1:end/2);
        nu = variables((end/2)+1:end);
        %% update
        pos = pos+uo;
        xhi = xhi + nu;
        if norm(G) <= 10^-10
            break;
        end
    end
        plot(pos(3:3:end),pos(1:3:end),'LineWidth',1,'color','m','LineStyle','-.');
        hold on;
%         plotter(xhi+xhi_ini,pos);

        xlim([-sum(len_ini),sum(len_ini)]);
        ylim([-sum(len_ini),sum(len_ini)]);
        drawnow;
        norm(G)
        iter
        
        pause(1);
%         clf;
end