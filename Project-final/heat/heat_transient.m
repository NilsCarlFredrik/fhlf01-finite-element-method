function [T_snap] = heat_transient(d0, K, C, F, snap_times, dt)
% [T_snaps]=heat_transient(d0, K, C, F, snap_times)
%-------------------------------------------------------------
% PURPOSE
%
%  Calculate the temperature distribution at given times using
%  the backward Euler scheme (step1 function from the calfem
%  toolbox).
%
% INPUT:
%         d0 :         Matrix [nnod x 1],     Initial temperature
%         K :          Matrix [nnod x nnod],  Stiffness matrix
%         C :          Matrix [nnod x nnod],  C-matrix
%         F :          Matrix[nnod x 1],      Force matrix
%         snap_times:  Matrix [1 x m],        Snap times
%                                             [time1 time2 ...]
%         dt :         Integer,               Time step size
%
% OUTPUT: T_snaps:      Matrix nnod x k
%-------------------------------------------------------------

%    Nils Broman, 2020-05-23
%-------------------------------------------------------------

%---------- Parameters used in step1 --------------
total_time = max(snap_times);        % Total time
theta = 1;                           % Specifies theta for backward euler
nsnap = size(snap_times,2);          % Number of snapshots
nhist = 0;                      % Number of history snaps saved (not used)
ip = [dt total_time theta nsnap nhist snap_times]; % Composed vector of
                                                   % inputs for step1

%----- Calculate T for each time step -------------
T_snap = step1(K,C,d0,ip,F); % Generate snapshots
