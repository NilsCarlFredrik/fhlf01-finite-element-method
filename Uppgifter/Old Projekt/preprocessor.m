function preprocessor(mesh)
% 
%
%
%
%
%

%------------ Initiate global constants -----------
global p e t coord nelm ndof dof edof
global th T_inf alpha_c T_0 Q Ex Ey
global E nu alpha rho cp k
%--------------------------------------------------

%------------------- Mesh -------------------------
%Loads mesh of points, edges and triangle elements (p,e,t) and places them
%in variables following the Calfem-format
load(mesh);
coord = p'/1000;            %Coordinates for nodal points [m]
t = t';                     %Triangle elements [node1 node2 node3 subdomain]
nelm = size(t,1);           %Number of elements
ndof = size(coord,1);       %Number of nodes
dof = (1:ndof);             %Nodal dof matrix
edof = [(1:nelm)' t(:,1:3)];%Topology matrix, elements with associated nodes

[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3); % Element node coordinates

%--------------------------------------------------

%--------------- General constants ----------------
th= 10e-3;               %Thickness
T_inf = 18;           %Outer uniform temperature
alpha_c = 40;         %
T_0 = 30;             %Initial temperature
Q = [0 0 0 5e7];   %Generated heat (only in Si)
%--------------------------------------------------

%--------------- Material constants ---------------
%Sets material constants and places them in global arrays formatted 
%after the subdomains in our mesh as follows: [Ag Cu Cu Si]

%Young's Modulus
E_Ag = 7;
E_Cu = 128;
E_Si = 165;
E = [E_Ag E_Cu E_Cu E_Si];

%Poisson's ratio
nu_Ag = 0.3;           
nu_Cu = 0.36;
nu_Si = 0.22;
nu = [nu_Ag nu_Cu nu_Cu nu_Si];

%Expansion coefficient
alpha_Ag = 4e-5;
alpha_Cu = 17.6e-6;
alpha_Si = 2.5e-6;
alpha = [alpha_Ag alpha_Cu alpha_Cu alpha_Si];

%Density
rho_Ag = 2500;
rho_Cu = 8930;
rho_Si = 2530;
rho = [rho_Ag rho_Cu rho_Cu rho_Si];

%Specific heat
cp_Ag = 1000;
cp_Cu = 386;
cp_Si = 703;
cp = [cp_Ag cp_Cu cp_Cu cp_Si];

%Thermal conductivity
k_Ag = 5;
k_Cu = 385;
k_Si = 149;
k = [k_Ag k_Cu k_Cu k_Si];
%--------------------------------------------------

