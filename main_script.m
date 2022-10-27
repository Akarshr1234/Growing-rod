input_scr;
%% Problem Setup
global L;               % gather operator ititialisation
L = zeros(npe*3,ndof,nel);                      % gather operator ititialisation

% gather operator definition
for el = 1:nel
    L(:,(el-1)*(3*(npe-1)) + 1 :((el-1)*(3*(npe-1))) + npe*3,el) = eye(npe*3);
end

% growth parameters
global Lo_dot;
global beta_growth; 
global gamma_growth; 
global dt;
global lg_growth;
global delta_B_stiffness;
global tau;

Lo_dot = 0.1;
beta_growth = 1; 
gamma_growth = 10/3;
lg_growth = 10;
tau = 1;

B_max_stiffness = 10*C(5,5);
Bo_stiffness = C(5,5);
delta_B_stiffness = B_max_stiffness - Bo_stiffness ;


len = len_finder(pos,nel);                  % length of each local
len_ini = len;
[gps,w] = gauss_legendre(ngp);              % gauss points and weights
    
dt = 1;
time = 50;
tsteps = time/dt;
timeins = linspace(0,time,tsteps);
% xlim([-sum(len_ini),sum(len_ini)]);
% ylim([-sum(len_ini),sum(len_ini)]);
% drawnow;
% pause(1);
len_iter = sum(len);
for t = 1:length(timeins)
    clc;
    fprintf('Time = %f \n',timeins(t))
    fprintf('Timestep = %d \n',t )
    pause(10);
    xhi = zeros(nnodes*3,1);
    len = len_finder(pos,nel);                  % length of each local
    variables = [pos;xhi];
%     variables = arc_length('residual_force',variables);
    variables = NR_iter('residual_force',variables);
    pos = variables(1:3*nnodes);
    xhi = variables(3*nnodes+1:end);
    [pos_ini,xhi_ini] = update_curvature(xhi_ini,pos_ini,xhi,pos); 
%     clf;
    plot(pos_ini(3:3:end),pos_ini(1:3:end),'r');
    hold on;
    plotter(xhi_ini,pos_ini,l/10,'r');
    hold on;
    plot(pos(3:3:end),pos(1:3:end),'b');
    xlim([-sum(len_ini),sum(len_ini)]*10);
    ylim([-sum(len_ini),sum(len_ini)]*10);
    plotter(xhi+xhi_ini,pos,l/10,'b');
    drawnow();
    pause(0.05);
    if t < length(timeins)
        clf;
    end
   

end
