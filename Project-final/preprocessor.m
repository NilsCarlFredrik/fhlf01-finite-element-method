function preprocessor(mesh)
% []=preprocessor(mesh)
%-------------------------------------------------------------
% PURPOSE
%
%  Initialize the parameters used within the program and stores
%  them in global variables to make them available to all other
%  functions.
%
% INPUT     mesh:        .mat,     Mesh used in the simulations
%-------------------------------------------------------------

%   Nils Broman, 2020-05-23
%-------------------------------------------------------------

%-------- Initialize global parameters -----------
global p e t coord enod nelm nnod dof dof_S edof edof_S Ex Ey % Mesh
global E nu alpha rho cp k Q_Si Q                             % Material
%-------------------------------------------------

%-------- Extract information from mesh ----------
load(mesh); % Load points, edges and triangles [p, e, t] from mesh
p = p/1000; % Transform point coordinates to meters

% Format information according to CALFEM standards
coord = p';             % Node coordinates
t = t';                 % Triangles, row format: [n1 n2 n3 sd]
                        %              (n=nodes, sd=subdomain)
enod = t(:,1:3);        % Nodes, row format: [n1 n2 n3]
e = e([1 2 5],:)';      % Reduce e and format rows: [n1 n2 edge_number]
nelm = size(enod,1);    % Number of elements
nnod = size(coord,1);   % Number of nodes
dof = (1:nnod)';        % Degrees of freedom for heat transfer,
                        % each row represents a node
dof_S = [dof*2-1 dof*2];% Dof for displacements

% Generate topology matrices adding element numbering
% to list of nodes associated with the elements
edof = [(1:nelm)' enod]; % For heat transfer

% Get element dof for x- and y-displacements
enod_y = enod*2;
enod_x = enod*2-1;

% Add enod_x and enod_y columnwise
edof_S = zeros(nelm, 7);    % Empty edof for displacements/stress
edof_S(:,1) = (1:nelm)';    % Fill first row with element numbers
edof_S(:,2:2:end) = enod_x; % x-axis dofs
edof_S(:,3:2:end) = enod_y; % y-axis dofs

[Ex,Ey]=coordxtr(edof,coord,(1:nnod)',3); % Element node coordinates
%-----------------------------------------

%--------------- Material constants ---------------
% Sets material constants and places them in global arrays formatted
% after the subdomain numbers in our mesh as follows: [Ag Cu Cu Si]

% Young's Modulus
E_Ag = 7;
E_Si = 165;
E_Cu = 128;
E = [E_Ag E_Cu E_Cu E_Si];

% Poisson's ratio
nu_Ag = 0.3;
nu_Si = 0.22;
nu_Cu = 0.36;
nu = [nu_Ag nu_Cu nu_Cu nu_Si];

% Expansion coefficient
alpha_Ag = 4e-5;
alpha_Si = 2.6e-6;
alpha_Cu = 17.6e-6;
alpha = [alpha_Ag alpha_Cu alpha_Cu alpha_Si];

% Density
rho_Ag = 2500;
rho_Si = 2530;
rho_Cu = 8930;
rho = [rho_Ag rho_Cu rho_Cu rho_Si];

% Specific heat
cp_Ag = 1000;
cp_Si = 703;
cp_Cu = 386;
cp = [cp_Ag cp_Cu cp_Cu cp_Si];

% Thermal conductivity
k_Ag = 5;
k_Si = 149;
k_Cu = 385;
k = [k_Ag k_Cu k_Cu k_Si];

% Heat generation
Q_Ag = 0;
Q_Si = 5e7;
Q_Cu = 0;
Q = [Q_Ag Q_Cu Q_Cu Q_Si];
%--------------------------------------------------
