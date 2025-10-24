function [sigma_eff_nodes] = stresses_stationary(Ed, T_stat)
% [sigma_eff_nodes]=stresses_stationary(Ed, T_stat)
%-------------------------------------------------------------
% PURPOSE
%
%  Calculate von Mises stresses on nodes due to displacements
%  caused by a temperature distribution.
%
% INPUT:  T_stat:            Matrix [nnod x 1], Temperature
%                                               destribution
%         Ed:                Matrix [nelm x 6], Element
%                                               displacements
%
% OUTPUT: sigma_eff_nodes:   Matrix [nnod x 1], Effective
%                                               nodal stress
%-------------------------------------------------------------

%   Nils Broman, 2020-05-25
%-------------------------------------------------------------

%------- Initialize global params --------
global th E nu alpha coord t nelm nnod T_0
%-----------------------------------------

sigma_eff_nodes = zeros(nelm,1);        % Empty vector of effectiv stress
sums_sigma_eff = zeros(nnod,2);         % Empty vector of sums, used to
                                        % calculate effectiv nodal stress

for elem=1:nelm
    sd = t(elem, 4);                % Element subdomain
    nodes = t(elem, 1:3);           % Element nodes
    ex = coord(nodes,1)';           % Node x-coordinates
    ey = coord(nodes,2)';           % Node y-coordinates
    D = hooke(2, E(sd), nu(sd));    % Constitutive matrix (2 for plane
                                    % strain, see hooke)


    es = plants(ex,ey,[2 th], D, Ed(elem,:)); % Calculate stresses
    % Compensate for thermal strains
    dT = mean(T_stat(nodes) - T_0);
    es = es - alpha(sd)*E(sd)*dT/(1-2*nu(sd))*[1 1 1 0];

    % Element Von Mises stress
    element_sigma_eff = sqrt(es(1)^2+es(2)^2+es(3)^2-es(1)*es(2) ...
                            -es(1)*es(3)-es(2)*es(3)+3*es(4)^2);

    % Increment first column, add von Mises stress to second (at rows
    % corresponding to element nodes).
    sums_sigma_eff(nodes,:) = sums_sigma_eff(nodes,:) + ...
                              [1 element_sigma_eff];

end

% Calculate nodal stresses as mean of surrounding elements
sigma_eff_nodes = sums_sigma_eff(:,2)./sums_sigma_eff(:,1);
