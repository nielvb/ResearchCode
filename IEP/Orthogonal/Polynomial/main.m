clear all
close all
addpath('../../') % to add givens and cgivens to path
m = 10;

% nodes
%z = [1:m]/20;
%z = cos([1:m]/m'*pi); %cos(rand(m,1)*pi); 
%z=rand(m,1) + 1i*rand(m,1);

%% equidistant on unit circle
int = linspace(0,2*pi,m+1);
int = int(1:m);
z = exp(1i*int);

%% weights
w = ones(m,1);

% Solve the Hessenberg IEP given above nodes and weights
[V,H] = Arnoldi(diag(z),w);
[Vup,Hup] = solveHIEP(z,w);
H = H(1:m,1:m);
V = V(:,1:m);
%% checks for solution to HIEP
norm(V(:,1)*norm(w)-w)
norm(diag(z)*V-V*H)

