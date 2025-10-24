function plot_temp(T, t_low, t_high, titles)
% []=plot_temp(T, t_low,t_high,titles)
%-------------------------------------------------------------
% PURPOSE
%  Plot one or more given temperature distributions.
%
% INPUT:  T :      Matrix nnod x k, where k is the number of 
%                  distributions to be plotted
%         t_low:   integer, lower bound for plot
%         t_high:  integer, upper bound for plot
%         titles:  1xk row vector with plot titles
%-------------------------------------------------------------

%------- Initialize global params --------
global Ex Ey edof
%-----------------------------------------

subplots = size(T,2);
for i = (1:subplots)
    % two subplots per row unless only 1 plot
    if (subplots == 1) 
        subplot(1,1,1);
    else
        subplot(round(subplots/2),2,i);
    end
    
    eT = extract(edof, T(:,i));
    fill(Ex', Ey', eT','EdgeColor','none');

    title(titles(i));
    
    % Set axes to battery dimensions
    axis equal
    ylim([-0.025 0.025])
    xlim([-0.0125 0.0125]);
    
    caxis([t_low t_high]);
    colormap(hot);
    colorbar;
    xlabel('x-position [m]');
    ylabel('y-position [m]');
end
