%% 1

%%%%%%%% Preprocessor %%%%%%%%%%%

A=10;
k=5;
Q=100;

a=2;
b=8;

Ta=0;
qb=15;

nen=2;
nelm=6;
ndof=nelm+1;

L=(b-a)/nelm;
D=A*k/L;

%%%%%%%%%%%%% Solver %%%%%%%%%%%%%%%%%%%

K=zeros(ndof);
f=zeros(ndof,1);
Ke=(A*k/L).*[1  -1; -1  1];

syms x
% Ne = 1/L*[-(1-xj)*x (1-xi)*x];
Be = 1/L*[1 -1];

for elnr=1:nelm
    Ne = 1/L*[-(x-(elnr+3)) (x-elnr)];
    fe = int(100*Ne',2+(elnr
    -1)*L,2+elnr*L);
    [K,f]=assem([elnr elnr elnr+1],K,Ke,f,fe);
end

fb = -qb*A*L;
f(end) = f(end)+fb
K
bc = [1 0];
a=solveq(K,f,bc)





















