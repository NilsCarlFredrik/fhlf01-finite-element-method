% Main script used to solve problems for the assignment in the course
% FHLF01 Finite Element Method 2020 at The Faculty of Engineering, Lund
% University
%
%  Nils Broman, 2020-05-26

%% Initialize
%------- Initialize global params --------------
clear
% Also variables not used in main for debugging purposes
global p e t coord enod nelm nnod dof dof_S edof edof_S Ex Ey % Mesh
global E nu alpha rho cp k Q_Si Q                             % Material
global alpha_c T_inf T_0 th                                   % Other
%-----------------------------------------------

%------------- General parameters --------------
alpha_c = 40;           % Heat transfer coefficient
T_inf = 18;             % Outer temperature [C]
T_0 = 30;               % Initial temperature [C]
th = 10e-3;             % Thickness [m]
%-----------------------------------------------

% Extract all information regarding grid and material parameters
preprocessor('mesh-refined')

%% a) 1 Solve stationary case for initial Q and reduced by 25%

% ---------------- Calculations  ----------
% Solve for initial Q
Q(4) = Q_Si;
[K_init,F_init] = heat_matrices();
T_init = solveq(K_init,F_init);

% Solve for reduced Q
Q(4) = Q_Si*0.75;
[K_red,F_red] = heat_matrices();
T_red = solveq(K_red,F_red);

% Max-/minimum temperatures for both cases
T_max_stat = max([T_init T_red])
T_min_stat = min([T_init T_red])
T_mean_stat = mean([T_init T_red])

%% a) 2 Plot stationary case for initial Q and reduced by 25%

%----------------- Plots -----------------------
% Rounded and scaled for visual purposes in plot
T_max_stat_plot = round(T_max_stat)*1.01;
T_min_stat_plot = round(T_min_stat)*0.99;

figure(2)

% Plots temperatures with initial Q (T_init), local scale
title_init = ['Stationary heat distribution for' newline ...
           'Q_S_i = 50 MW/m^3'];
subplot(2,2,1)
plot_temp(T_init, T_min_stat_plot(1), T_max_stat_plot(1), title_init);

% Plots temperatures with reduced Q (T_red), local scale
title_red = string(['Stationary heat distribution for' newline ...
           'Q_S_i = 37.5 MW/m^3']);
subplot(2,2,2)
plot_temp(T_red, T_min_stat_plot(2), T_max_stat_plot(2), title_red);

% Plots temperatures with initial Q (T_init), global scale
subplot(2,2,3)
plot_temp(T_init, T_min_stat_plot(2), T_max_stat_plot(1), '');

% Plots temperatures with reduced Q (T_red), global scale
subplot(2,2,4)
plot_temp(T_red, T_min_stat_plot(2), T_max_stat_plot(1), '');


%% b) 1 Solve transient solutions for initial and reduced Q

% -------- General transient parameters --------
T_start = T_0*ones(nnod,1);         % Starting temperature
snap_times = 60.*[5 15 20];% Snapshot times (minutes translated to seconds)
dt = 1;                         % Step size (in seconds)

%---------- Calculations initial Q -------------
Q(4) = Q_Si; % Reset to initial Q
[K_init,F_init,C_init] = heat_matrices();

% Transient solutions at times given by snap_times
T_trans_init = heat_transient(T_start, K_init, C_init, F_init, ...
                              snap_times, dt);
T_all_init = [T_start T_trans_init];

% Max-/minimum temperatures for each time (rounded to 1 decimal)
T_max_init = round(max(T_all_init),1)
T_min_init = round(min(T_all_init),1)
%-----------------------------------------------

%-------------- Calculations reduced Q ---------
Q(4) = Q_Si*0.75; % Reduce used Q by 25%
[K_red,F_red,C_red] = heat_matrices();

% Transient solutions at given times
T_trans_red = heat_transient(T_start, K_red, C_red, F_red, snap_times, dt);
T_all_red = [T_start T_trans_red];

% Max-/minimum temperatures for each time (rounded to 1 decimal)
T_max_red = round(max(T_all_red),2)
T_min_red = round(min(T_all_red),20)
%-----------------------------------------------

%% b) 2 Transient plots

% Scaled max/min for all times for visual purposes in plot
T_max_trans_plot = max([T_max_init T_max_red])+0.5;
T_min_trans_plot = min([T_min_init T_min_red])-0.5;

%------------- Plots init trans ----------------
figure(3)
titles = ['    Start: T_m_a_x = ' + string(T_max_init(1)) + ' C';
          '    5 min: T_m_a_x = ' + string(T_max_init(2)) + ' C';
          '    12 min: T_m_a_x = ' + string(T_max_init(3)) + ' C';
          '    20 min: T_m_a_x = ' + string(T_max_init(4)) + ' C'];

for(i=1:4)
    subplot(4,1,i)
    plot_temp(T_all_init(:,i), T_min_trans_plot, T_max_trans_plot, ...
              titles(i));
end
sgtitle('Q_S_i = 50 MW/m^3', 'FontSize', 12)

%------------ Plots red trans ------------------
figure(4)
titles = ['    Start: T_m_a_x = ' + string(T_max_red(1)) + ' C';
          '    5 min: T_m_a_x = ' + string(T_max_red(2)) + ' C';
          '    12 min: T_m_a_x = ' + string(T_max_red(3)) + ' C';
          '    20 min: T_m_a_x = ' + string(T_max_red(4)) + ' C'];

for(i=1:4)
    subplot(4,1,i)
    plot_temp(T_all_red(:,i), T_min_trans_plot, T_max_trans_plot, ...
              titles(i));
end
sgtitle('Q_S_i = 37.5 MW/m^3', 'FontSize', 12)

%% b) Time till stationary
% I realized this was part of the assignment too late hence the extremely
% inelegant approach.

% -------- General transient parameters --------
sim_time= 52*24*60;                   % Simulated time (in m?nutes)
snap_times = 60*(1:sim_time);         % Snapshot times (one every minute)
dt = 10;                              % Step size (10 seconds)

% ------------ Initial Q -----------------------

% Transient solutions at times given by snap_times
T_trans_init = heat_transient(T_start, K_init, C_init, F_init, ...
                              snap_times, dt);

for(i=1:sim_time)
    if (T_init-T_trans_init(:,i)) < 0.1*ones(length(T_init),1) ...
        == ones(length(T_init),1);
        i/60
        break
    end
end

% ------------ Reduced Q -----------------------

% Transient solutions at given times
T_trans_red = heat_transient(T_start, K_red, C_red, F_red, snap_times, dt);

for(i=1:sim_time)
    if (T_init-T_trans_red(:,i)) < 0.1*ones(length(T_init),1) ...
        == ones(length(T_init),1);
        i/60
        break
    end
end

'Finished'


%% c) 1 Calculate displacements and Von Mises stress for initial Q
%%%%%%%%% RUN a) 1 BEFORE THIS, needs T_init %%%%%%%%%%%%

% Calculates diplacement matrices
[K_init, F0_init, bc] = disp_matrices(T_init);

u_init = solveq(K_init,F0_init,bc); % Solves displacements [x1,y1,x2,y2...]
u_x_init = u_init(1:2:end,:);       % Extracts x-displacements
u_y_init = u_init(2:2:end,:);       % Extracts y-displacements
max_disp_init = max(sqrt(u_x_init.^2+u_y_init.^2)); % Max node displacement

Ed_init = extract(edof_S,u_init);

sigma_eff_nodes_init = stresses_stationary(Ed_init, T_init);
[max_stress_init,max_stress_node_init] = max(sigma_eff_nodes_init);



%% c) 2 Calculate displacements and Von Mises stress for reduced Q
%%%%%%%%% RUN a) 1 BEFORE THIS, needs T_red %%%%%%%%%%%%

% Calculates diplacement matrices
[K_red, F0_red, bc] = disp_matrices(T_red);

u_red = solveq(K_red,F0_red,bc);    % Solves displacements [x1,y1,x2,y2...]
u_x_red = u_red(1:2:end,:);         % Extracts x-displacements
u_y_red = u_red(2:2:end,:);         % Extracts y-displacements
max_disp_red = max(sqrt(u_x_red.^2 + u_y_red.^2)); % Max node displacement

Ed_red = extract(edof_S,u_red);

sigma_eff_nodes_red = stresses_stationary(Ed_red, T_red);
[max_stress_red,max_stress_node_red] = max(sigma_eff_nodes_red);

%% c) 3 Plots displacements and stresses for both initial and reduced Q

titles = [string(['Displacement and Von Mises stress' newline ...
                 'for Q_S_i = 50 MW/m^3' newline]);
          string(['Displacement and Von Mises stress' newline ...
                 'for Q_S_i = 37.5 MW/m^3' newline])];

% Plot for initial Q
% subplot(2,1,1);
figure(10)
plot_stress_disp(u_init, sigma_eff_nodes_init, max_stress_node_init, ...
                 titles(1))
% Plot for reduced Q
% subplot(2,1,2);
figure(11)
plot_stress_disp(u_red, sigma_eff_nodes_red, max_stress_node_red, ...
                titles(2))
