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
global age_vector;   % age of an element out of growing region

global pos_ini;      % undeformed configuration position
global xhi_ini;      % undeformed configuration rotation
global Force;        % natural boundary condition
global body_Force;   % body force vector
global len_ini;      % initial length vector


l = 5;                         % total length
nel = 20;                       % number of locals
npe = 2;                        % number of nodes per local
nnodes = (npe-1)*nel + 1;       % total number of nodes
ndof = 3* nnodes;               % total number of positional degree of freedom
ngp = 3;                        % number of gauss points to be used

% n_bar and m_bar value
n_bar_c = [-1; 0; 0];            % constant body force vector
m_bar_c = [0; 0; 0];            % constant body moment vector



n_bar = zeros(nnodes*3,1);      % initialisation of body force
m_bar = zeros(nnodes*3,1);      % initialisation of body moment



%% Essential Boundary Conditions
bcs = [1,2,3,nnodes*3+1,nnodes*3+2,nnodes*3+3];%,6*nnodes,6*nnodes-2,6*nnodes-1 ,nnodes*3-2,nnodes*3-1,nnodes*3]; %];%  bounadary condition dof
bc_val = [0,0,0,0,0,0,0,0,0,0,0,0]; % boundary condition values

% body forces and body moment definition
for node = 1:nnodes
    n_bar(((node-1)*3+1):((node-1)*3+3)) = n_bar_c;
    m_bar(((node-1)*3+1):((node-1)*3+3)) = m_bar_c;
end




% C = diag([10^7 10^8 10^6 10^4 10^5 10^3]);
% C = diag((1/pi^2)*ones(6,1));                        % material parameters
C = diag((100/8)^3*ones(6,1));     
C(1,1) = 100*((100/8)^3);
C(2,2) = C(1,1);
C(3,3) = C(1,1);
Be_stiffness = C(5,5)*ones(nel,1);
% C = diag([1.61538 *10^8,1.61538 *10^8,3.5 *10^7,3.5 *10^7,3.5 *10^7,3.5 *10^7]);
% C = diag([3600*10^10,3600*10^10,7200*10^10,7200*0.5,7200*0.5,7200]);
age_vector = zeros(nel,1);
% initial position discription
% posz = linspace(0,l,nnodes)';
posz = mesher(0.75,nel,l);
pos = zeros(nnodes*3,1);
for node = 1:nnodes
   pos(3*(node)) = posz(node); 
end
% pos = pos_arc*10;
% pos(1:3:end) = -pos(1:3:end);
pos_ini = pos;
% initial rotation description

xhi = zeros(nnodes*3,1);
% xhi_ini = xhi_arc;
xhi_ini = zeros(nnodes*3,1);
len_ini = len_finder(pos_ini,nel);

% 
% posz = linspace(0,1,nnodes)';
% posy = 0.0*(-1-cos(2*pi*(linspace(0,1,nnodes))))';
% % posz = mesh(0.9,nel,l);
% pos = zeros(nnodes*3,1);
% for node = 1:nnodes
%    pos(3*(node)) = posz(node); 
% end


%% imposing natural boundary condition
Force = zeros(nnodes*6,1);
% Force(((nnodes+1)/2)*3,1) = -400;
% Force(((nnodes+1)/2)*3-2,1) = 0.2;
% Force(nnodes*6-1,1) = 2*pi*C(4,4)/l;



%% body force
body_Force = zeros(nnodes*6,1);
