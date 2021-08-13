function [Qtil, Htil,counter]= mergeIEP_HH(Q1,Q2,H1,H2,w1,w2)
% Merges two solutions to a Hessenberg IEP to form solution to a larger Hess IEP using Householder reflectors
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
for i=1:n-1
  irow = min(i+1,m); % in case k<n, then the length of the vector to eleminate is upper bounded
  x = [Hhat(i+1,i);Hhat([n+1:n+irow],i)];
  omega = angle(x(1));
  e = exp(1i*omega);
  v = (x-e*norm(x)*[1;zeros(irow,1)])/(norm(x-e*norm(x)*[1;zeros(irow,1)]));
  Pdot = eye(irow+1)-2*v*v';
  P = eye(n+m);
  P([i+1,n+1:n+irow],i+1) = Pdot(:,1);
  P([i+1,n+1:n+irow],[n+1:n+irow]) = Pdot(:,2:irow+1);
  Qhat = Qhat*P';
  Hhat = P*Hhat*P';
  Hhat([n+1:n+irow],i) = 0; % set explicitly to zero
end

%% Restore Hessenberg structure in trailing mxm principa submatrix
%for i=1:m-1
%  for j = m:-1:i+1
%    P= I; P([n+j-1,n+j],[n+j-1,n+j]) = givens(Hhat(n+j-1,n+i-1),Hhat(n+j,n+i-1))';
%    counter = counter + 1;%%%%%%%%%%%%%%%
%    
%    Qhat = Qhat*P;
%    Hhat = P'*Hhat*P;
%    Hhat(n+j,n+i-1) = 0; % set explicitly to zero
%  end
%end
%ifinal = min(n+n+1,n+m);
%for i=n:ifinal-1
%i
%  x = Hhat([i+1:ifinal],i);
%  omega = arg(x(1));
%  e = exp(1i*omega);
%  v = (x-e*norm(x)*[1;zeros(ifinal-i-1,1)])/(norm(x-e*norm(x)*[1;zeros(ifinal-i-1,1)]));
%  Pdot = eye(ifinal-i)-2*v*v';
%  P = eye(n+m);
%  P([i+1:ifinal],[i+1:ifinal]) = Pdot;
%  Qhat = Qhat*P';
%  Hhat = P*Hhat*P';
%  %Hhat([n+1:n+irow],i) = 0; % set explicitly to zero
%end


ifinal = n+m;
i = n;
while i<n+m-1
  if n+i+1<=n+m
    ifinal = n+i+1;
  end 
  x = Hhat([i+1:ifinal],i);
  omega = angle(x(1));
  e = exp(1i*omega);
  v = (x-e*norm(x)*[1;zeros(ifinal-i-1,1)])/(norm(x-e*norm(x)*[1;zeros(ifinal-i-1,1)]));
  Pdot = eye(ifinal-i)-2*v*v';
  P = eye(n+m);
  P([i+1:ifinal],[i+1:ifinal]) = Pdot;
  Qhat = Qhat*P';
  Hhat = P*Hhat*P';
  Hhat([i+2:ifinal],i) = 0; % set explicitly to zero
  
  i = i + 1;
end


Qtil = Qhat;
Htil = Hhat;


return;
