%% 3.4 Ex

clear all
close all

% addpath C:\Users\Nils\Google Drive\Skolarbete\Lth\Ar 3\FHLF01 Finita Element Metoden\Matlab\calfem

%%%%%%%%%%% Preprocessor
% Definerar parametrar

A=1; 
E=1;
L=1;

% Jag v?ljer att l?gga nod med frihetsgraderna u_5 u_6 i origo

Coord = [ -L            0
          -L/sqrt(2)    L/sqrt(2)
          0             0
          L/sqrt(2)     L/sqrt(2)   ];
      
Edof = [1       1 2 5 6
        2       3 4 5 6
        3       5 6 7 8];
    
nen = 2;    % # noder per element
ndof = 8;   % # frihetsgrader
nelm = 3;   % # elements

Dof= [ 1 2
       3 4
       5 6
       7 8];
[Ex, Ey] = coordxtr(Edof, Coord, Dof, nen);

%%%%%%%%%%%%%%%%%%%%%%%% Solver

K=zeros(ndof); % styvhetsstorlek som frihetsgrader

for elnr=1:nelm
    % Brukar vara n?got h?r som ?r specifikt f?r elementet E = minaE(elnr);
    % (typ materialegenskaper som tr?, plast osv.)
    Ke=bar2e(Ex(elnr,:),Ey(elnr,:),[A E]);
    K=assem(Edof(elnr,:),K,Ke);    
end

%% 3.4 Ex ut?kad (mer ?n uppgiften)

clear all
close all

%addpath /...

%%%%%%%%%%% Preprocessor
% Definerar parametrar

A=1; 
E=1;
L=1;
Fl=10; %%%%%%%


% Jag v?ljer att l?gga nod med frihetsgraderna u_5 u_6 i origo

Coord = [ -L            0
          -L/sqrt(2)    L/sqrt(2)
          0             0
          L/sqrt(2)     L/sqrt(2)   ];
      
Edof = [1       1 2 5 6
        2       3 4 5 6
        3       5 6 7 8];
    
nen = 2;    % # antal noder per element
ndof = 8;   % # frihetsgrader
nelm = 3;   % # elements

Dof= [ 1 2
       3 4
       5 6
       7 8];

F = zeros(ndof,1); %%%%%%

F(5)=Fl; %%%%%kraft i frihetsgrad 5
F(6)=-Fl; %%%%%%

[Ex, Ey] = coordxtr(Edof, Coord, Dof, nen);

%%%%%%%%%%%%%%%%%%%%%%%% Solver

K=zeros(ndof); % styvhetsmatris storlek som frihetsgrader

for elnr=1:nelm
    % Brukar vara n?got h?r som ?r specifikt f?r elementet E = minaE(elnr);
    % (typ materialegenskaper som tr?, plast osv.)
    Ke=bar2e(Ex(elnr,:),Ey(elnr,:),[A E]);
    K=assem(Edof(elnr,:),K,Ke);    
end


bc=[ 1 0
     2 0
     3 0
     4 0
     7 0
     8 0 ];  %l?s fast frihetsgrader (insp?nningar)
 
 a=solveq(K,F,bc); %f?rskjutningar nod ?
 
 %%%%%%%%%%%%%%%%%%%%%%%%%% Extra
 
 
 ed=extract(Edof,a); %f?rskjutningar f?r element
 
 es=bar2s(Ex,Ey,[A E],ed); %sp?nningar
 
 plotpar=[1 4 1]; %hur linjen ska se ut
 eldraw2(Ex,Ey,plotpar)
 
 plotpar=[1 2 1]; %hur linjen ska se ut
 eldisp2(Ex,Ey,ed,plotpar)
 
 
 
 
 
 
 
 
 
 