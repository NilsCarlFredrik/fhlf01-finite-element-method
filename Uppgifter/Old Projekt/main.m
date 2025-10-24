%% Set parameters

clear
%------------ Initiate global constants -----------
global p e t coord nelm ndof edof Ex Ey    %Element parameters
global th T_inf alpha_c T_0 Q              %General constants
global E nu alpha rho cp k                 %Material constants

preprocessor('testmesh2.mat')

[K, F, C] = heat_matrices();

T_stat = solve(K,F);

eT = extract(edof, T_stat);
 
% length(coord(:,1))
% length(coord(:,2))
% length(eT)

fill(Ex', Ey', eT');

%%
T_low = 10;
T_high = 20;

title = string(['Stationary heat distribution for' newline ...
           'T_i_n_f = 15 ?C, Q_c_o_r_e = 10^5 W/m^3'])

plot_temp(T_stat, T_low, T_high, title);
peak_temp_stat = max(T_stat)