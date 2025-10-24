function [K, F, C] = heat_matrices()
%
%
%
%
%
%

%---- Used global constants (see preprocessor) ----
global ndof nelm t coord k Q th edof rho cp e T_inf alpha_c

%------------ Create empty matrices ---------------
K = zeros(ndof);            %Empty stiffness matrix
fl = zeros(ndof, 1);        %Empty load force vector
fb = zeros(ndof, 1);        %Empty boundary force vector
fc = zeros(ndof);
C = zeros(ndof);            %Empty 

%-------------- Compute K, fl & C -----------------
for elnr = (1:nelm)
   elnod = t(elnr, 1:3);                 %Element nodes
   sd = t(elnr,4);                       %Element sub domain
   ex = coord(elnod, 1)';                %Element node x coordinates
   ey = coord(elnod, 2)';                %Element node y coordinates
   D = k(sd)*eye(2);                     %Constitutive matrix
   Qe = Q(sd);                           %Generetad heat in element
   
   [Ke,fle] = flw2te(ex,ey,th,Qe);       %Compute local K & fl for element
%    [K,fl]=assem(edof(elnr,:),K,Ke,fl,fle);       %Assemble K & fl
   
   Ce = plantml(ex,ey,rho(sd)*cp(sd));   %Compute element C
%    C = assem(edof(elnr,:),C,Ce);                 %Assemble C

%%%%%%%%%%%% Jonas assemblering %%%%%%%%%%%
%     indx = edof(elnr, 2:end);
%     K(indx, indx) = K(indx, indx) + Ke;
%     C(indx, indx) = C(indx, indx) + Ce;
%     fl(indx) = fl(indx) + fle;
end
size(Ke)
size(fle)
%----------- Compute convection on boundaries --------
% conv_bound = [1 2 9];               %Boundaries with convection
conv_bound = [2 3 9];               %Boundaries with convection (2)
i = ismember(e(5,:), conv_bound);   %Logical array of edges with convection
e_conv = e(1:2,i);                  


for enr = (1:length(e_conv))
    edge = e_conv(:,enr);               %Edge nodes
    ex = coord(edge,1);                 %x coordinates for edge nodes
    ey = coord(edge,2);                 %y coordinates for edge nodes
    L = sqrt(diff(ex)^2 + diff(ey)^2);  %Edge length
    
    fbce = th*L*alpha_c*T_inf/2*[1;1];  %Element convection force boundary vector
    fkce = th*L*alpha_c/6*[2 1; 1 2];   %Element convection boundary to K
    
%     n1=edge(1)                         %Node 1
%     n2=edge(2)                         %Node 2
%     'asdf'
%     fb([n1,n2]) = fb([n1,n2]) + fbce;
%     K([n1 n2], [n1 n2]) = K([n1 n2], [n1 n2]) + fkce;

    fbce = th*L*alpha_c*T_inf/2*[1;1;0];  %Element convection force boundary vector
    fkce = th*L*alpha_c/6*[2 1 0;
                           1 2 0;
                           0 0 0]          %Element convection boundary to K
    [fc, fb] = assem(edof,fc,fkce,fb,fbce);
end

% K = K + fc;
F = fl + fb;



















