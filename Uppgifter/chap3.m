%% 3.1

%%%%%%%%%%%% pre processor

k = 8;

nen = 2;
ndof = 4;
nelm = 5;

Edof = [1  1 2
        2  1 3
        3  2 4
        4  2 3
        5  3 4];
    
Ke = spring1e(k);

K = zeros(ndof);
F = zeros(ndof,1);

F(2) = 0;
F(3) = 20;

bc = [1 1
      4 0];
  
%%%%%%%%% Solver

for elnr=1:nelm
    K = assem(Edof(elnr,:),K,Ke);
end

u = solveq(K,F,bc);

%%%%%%%%%% Post processor (funkar inte, kolla upp)

ed=extract(Edof,u);
es=spring1s(k,ed)

%% 3.2

%%%%%%%%%%%% Preprocessor

k = 10;

nen=2;
ndof=3;
nelm=4;

%Dof = [ 1
 %       2
  %      3];

Edof = [1  1 2
        2  1 2
        3  2 3
        4  1 3];

%%%%%%%%%%%%% solver

K = zeros(ndof);
F = zeros(ndof,1);


for elnr=1:nelm
    Ke = k.*[1 -1; -1 1];
    K=assem(Edof(elnr,:),K,Ke);
end

F(3) = 10;

u = solveq(K,F,[1 0]);

%% 3.3 

Ke = [1 -1; -1 1];

nen = 2;
ndof = 4;
nelm = 5;

Dof = [];

Edof = [1  1 2
        2  1 3
        3  2 3 
        4  3 4
        5  2 4];

K = zeros(ndof);
I = zeros(ndof,1);

for elnr=1:nelm
    K=assem(Edof(elnr,:),K,Ke);
end

bc = [2 0
      4 0];
  
I(1) = 1;


V = solveq(K,I,bc);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  