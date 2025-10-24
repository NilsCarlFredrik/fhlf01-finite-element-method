function [K, F0, bc] = disp_matrices(T_stat)
% [K, F0, bc]=disp_matrices(T_stat)
%-------------------------------------------------------------------
% PURPOSE
%
%   Calculates displacement matrices for a system with given
%   temperature distribution, uses global variables.
%
% INPUT:  T_stat :      Matrix nnod x 1,           Temperature
%                                                  distribution vector
%
% OUTPUT: K:            Matrix [2*nnod x 2*nnod],  Stiffness matrix
%         F0:           Matrix [2*nnod x 1],       Temperature strain
%                                                  vector
%         bc:           Matrix [locked nodes x 2], Boundary condition
%                                                  vector
%------------------------------------------------------------------

%   Nils Broman, 2020-05-24
%------------------------------------------------------------------

%--- Initialize global parameters --------
global th E nu alpha coord t e nelm nnod edof_S T_0
%-----------------------------------------

%------- Generate K & F0 -----------------
K = zeros(2*nnod);                  % Empty stiffness matrix
Fb = zeros(2*nnod,1);               % Empty boundary vector
F0 = zeros(2*nnod,1);               % Empty temperature strain vector

for elem = (1:nelm)                 % Iterate over all elements
    sd = t(elem, 4);                % Element subdomain
    nodes = t(elem, 1:3);           % Element nodes
    ex = coord(nodes,1)';           % Node x-coordinates
    ey = coord(nodes,2)';           % Node y-coordinates

    D = hooke(2, E(sd), nu(sd));    % 2 => Plane strain
    D = D([1 2 4],[1 2 4]);         % Remove terms associated with sigma_zz

    % Generate element stiffness matrix
    Ke = plante(ex, ey,[2 th], D);

    %Calculate thermal strains contribution
    dT = mean(T_stat(nodes)) - T_0; % T of element taken as mean of nodes
    [Be, A] = planteBe(ex, ey);     % Generate B matrix and element area
    ep_dT = (1+nu(sd))*alpha(sd)*dT*[1;1;0];
    f0e = Be'*D*ep_dT*A;            % Element thermal load vector

    % Assemble element force vector and stiffness
    % matrix into global K and Fl
    indx = edof_S(elem, 2:end);
    K(indx, indx) = K(indx, indx) + Ke;
    F0(indx) = F0(indx) + f0e;
end
%-----------------------------------------

%------- Generate bc ---------------------
% Edges with boundary conditions (from mesh)
fixed_segments_x = [4 20 21 22 23];
fixed_segments_y = [4];

% Add conditions for fixed x-dofs
indx = ismember(e(:,3), fixed_segments_x);
edge_nodes_fixed = unique([e(indx,1); e(indx,2)]);
bc_x = [ 2*edge_nodes_fixed-1    zeros(size(edge_nodes_fixed))  ];
% 2*x-1 compensates x-y-stacking

% Add conditions for fixed y-dofs
indx = ismember(e(:,3),fixed_segments_y);
edge_nodes_fixed = unique([e(indx,1); e(indx,2)]);
bc_y = [ 2*edge_nodes_fixed    zeros(size(edge_nodes_fixed))  ];
% 2*x compensetes x-y-stacking

bc = sort([bc_x ; bc_y]);
%-----------------------------------------
