function [K, F, C] = heat_matrices()
% [K, F, C]=heat_matrices()
%-------------------------------------------------------------
% PURPOSE
%
%  Calculate temperature distribution matrices for stationary
%  and transient problems. Uses global variables, hence no input.
%
% OUTPUT: K :      Matrix nnod x nnod
%         F :      Matrix nnod x 1
%         C :      Matrix nnod x nnod
%-------------------------------------------------------------

%    Nils Broman, 2020-05-23
%-------------------------------------------------------------

%------------ Initialize global params ----------
global alpha_c T_inf th k cp rho Q coord t e nnod edof

%------------ Calculate K, Fl & C ---------------
K = zeros(nnod);                % Empty stiffness matrix
Fl = zeros(nnod,1);             % Empty load vector
Fb = zeros(nnod,1);             % Empty boundary vector
Kc = zeros(nnod);               % Empty convection matrix
C = zeros(nnod);                % Empty C matrix

for elnr = (1:size(t,1))        % Iterate over all elements
    sd = t(elnr, 4);            % Element subdomain
    nodes = t(elnr, 1:3);       % Element nodes
    ex = coord(nodes,1)';       % X coordinates of element nodes
    ey = coord(nodes,2)';       % Y coordinates of element nodes
    D = eye(2)*k(sd);           % Constitutive matrix

    % Calculate element stiffness matrix and load vector
    [Ke, fle] = flw2te(ex, ey, th, D, Q(sd));
    % Assemble Ke and fle in global matrices
    [K,Fl]=assem(edof(elnr,:),K,Ke,Fl,fle);

%     ce = cp(sd);
%     rhoe = rho(sd);
    Ce = plantml(ex,ey,rho(sd)*cp(sd));
    C = assem(edof(elnr,:),C,Ce);                 % Assemble C

end
%-----------------------------------------------

%------- Calculate boundary forces -------------
conv_boundary = [2 3 9];                % Boundaries with convection
indx = ismember(e(:,3),conv_boundary); % Logical array of edges on boundary
edges_conv = e(indx,:);      % All edges belonging to convection boundaries

for edge_i = (1:length(edges_conv))  % Iterate ofer edges on boundary
    edge = edges_conv(edge_i,:);

    n1 = edge(1);             % 1st node (on boundary)
    n2 = edge(2);             % 2nd node (on boundary)

    ex = coord(edge(1:2),1);  % x coordinates for edge nodes
    ey = coord(edge(1:2),2);  % y coordinates for edge nodes
    dx = diff(ex);
    dy = diff(ey);

    L = sqrt(dx^2+dy^2);            % Edge length
    fbe = L*alpha_c*T_inf/2*th*[1;1];  % Element boundary force for F
    Kce = L*alpha_c/6*th*[2 1; 1 2];   % Element boundary for K

    % Assemble Kc and Fb
    Fb([n1,n2]) = Fb([n1,n2]) + fbe;                    % Assemble Fb
    Kc([n1 n2], [n1 n2]) = Kc([n1 n2], [n1 n2]) + Kce;  % Assemble new K
                                                        % with convection

end
%------------------------------------------------

F = Fl + Fb;        % Add Fl and Fb for total F
K = K + Kc;         % Add convection to K
