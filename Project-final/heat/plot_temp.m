function plot_temp(T, T_low, T_high, the_title)
% []=plot_temp(T, T_low, T_high, titles)
%-------------------------------------------------------------
% PURPOSE
%
%  Plots temperature distribution, uses global variables.
%
% INPUT:  T :           Matrix [nnod x 1], Temperature vector
%         T_low:        Integer,           Lower bound for plot
%         T_high:       Integer,           Upper bound for plot
%         the_title:    String,            Plot title
%-------------------------------------------------------------

%  Nils Broman, 2020-05-23
%-------------------------------------------------------------


%---------- Global Parameters ------------
global Ex Ey edof
%-----------------------------------------

eT = extract(edof,T);

% Extend arrays to mirror plot in y-axis
Ex_plot = [-flip(Ex') Ex'];
Ey_plot = [flip(Ey') Ey'];
eT_plot = [flip(eT') eT'];

fill(Ex_plot, Ey_plot, eT_plot, 'EdgeColor', 'white');

% title(the_title, 'FontSize', 8)
title(the_title);

axis equal
ylim([-0.0005 0.0075])
xlim([-0.0045 0.0045]);

caxis([T_low T_high]);
colormap(hot);
colorbar;
xlabel('x-position [m]');
ylabel('y-position [m]');
