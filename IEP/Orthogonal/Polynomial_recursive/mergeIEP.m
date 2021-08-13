function [Qtil, Htil,counter]= mergeIEP(Q1,Q2,H1,H2,w1,w2)
% Merges two solutions to a Hessenberg IEP to form solution to a larger Hess IEP
% INPUT:
%       Q1,Q2 = unitary matrix of solution of some Hess IEP
%       H1,H2 = Hessenberg matrix of solution of some Hess IEP
%       w1,w2 = weigths of the Hess IEPs, these are required because the scaling in Q1 and Q2 is unknown
% OUTPUT:
%       Qtil,Htil = solution to Hess IEP obtained by combination of two given IEPs
% Note: dimensions of problems do not need to match

n = length(Q1);
m = length(Q2);
I = eye(n+m);

% Merge the given matrices and weigts
Qhat = [Q1, zeros(n,m);...
        zeros(m,n),Q2];
Hhat = [H1, zeros(n,m);...
        zeros(m,n),H2];
wtil = [w1;w2]; 
   
% Enforce orthogonality condition wrt new weigth vector
P= I; P([1,n+1],[1,n+1]) = givens(norm(w1),norm(w2))';
counter = 1;%%%%%%%%%%%%%%%

Qhat = Qhat*P;
Hhat = P'*Hhat*P;

% Enforce structure on Hhat
%% annihilate nonzeros in first n columns
for i=2:n
  for j=min(m,i):-1:2
    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Hhat(n+j-1,i-1),Hhat(n+j,i-1))';
    counter = counter + 1;%%%%%%%%%%%%%%%
    
    Qhat = Qhat*P;
    Hhat = P'*Hhat*P;  
    Hhat(n+j,i-1) = 0; % set explicitly to zero
  end
  P= I; P([i,n+1],[i,n+1]) = givens(Hhat(i,i-1),Hhat(n+1,i-1))';
  counter = counter + 1;%%%%%%%%%%%%%%%
  
  Qhat = Qhat*P;
  Hhat = P'*Hhat*P;
  
  Hhat(n+1,i-1) = 0; % set explicitly to zero
end

%% Restore Hessenberg structure in trailing mxm principa submatrix
for i=1:m-1
  for j = m:-1:i+1
    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Hhat(n+j-1,n+i-1),Hhat(n+j,n+i-1))';
    counter = counter + 1;%%%%%%%%%%%%%%%
    
    Qhat = Qhat*P;
    Hhat = P'*Hhat*P;
    Hhat(n+j,n+i-1) = 0; % set explicitly to zero
  end
end

Qtil = Qhat;
Htil = Hhat;


return;
