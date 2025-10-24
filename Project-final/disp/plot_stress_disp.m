function plot_stress_disp(u, sigma_eff_nodes, max_stress_node,the_title)
% []=plot_stress_disp(u, sigma_eff_nodes, max_stess_node, the_title)
%-------------------------------------------------------------
% PURPOSE
%
% Plots displacements and stresses for u, uses global variables.
%
% INPUT:  u:                     Matrix [2*nnod x 1],   Displacements
%         sigma_eff_nodes :      Matrix [nnod x 1],     Node stresses
%         max_stress_node:       Integer,               Node number of
%                                                       maximum stress
%         the_title:             String,                Title of plot
%-------------------------------------------------------------

%   Nils Broman, 2020-05-24
%-------------------------------------------------------------

%------- Initialize global params --------
global Ex Ey edof coord edof_S
%-----------------------------------------

% Calculate displacements of elements
Ed = extract(edof_S,u);

% Add displacements to nodes
Exd = Ex + Ed(:,1:2:end);
Eyd = Ey + Ed(:,2:2:end);
% Extracts sress on nodes
Esigma = extract(edof, sigma_eff_nodes);

% Extends arrays for plot to mirror solution around y-axis
Exd_plot = [-flip(Exd) ; Exd];
Eyd_plot = [flip(Eyd) ; Eyd];
Esigma_plot = [flip(Esigma) ; Esigma];

hold off
% Plots new nodes and element stresses
fill(Exd_plot',Eyd_plot',Esigma_plot','EdgeColor','white', ...
    'HandleVisibility','off')
hold on

% Calculate coordinates for node of maximum stress
max_stress_node_x = coord(max_stress_node,1) + u(max_stress_node*2-1);
max_stress_node_y = coord(max_stress_node,2) + u(max_stress_node*2);

% Plots points on nodes of maximum stress
fill_color = 'black';
edge_color = fill_color;

plot(max_stress_node_x*[-1 1],max_stress_node_y*[1 1],'o','MarkerSize', ...
    3,'MarkerFaceColor',fill_color, 'MarkerEdgeColor',edge_color);
plot(max_stress_node_x*[-1 1],max_stress_node_y*[1 1],'o','MarkerSize', ...
    10,'MarkerEdgeColor',edge_color);

legend(['Points of maximum stress'])

% Sets color of stress plot and adds colorbar to plot
caxis([0 13])
colormap(jet);
colorbar;
axis equal
ylim([-0.0005 0.0085])
xlim([-0.005 0.005])
title(the_title)

hold off
