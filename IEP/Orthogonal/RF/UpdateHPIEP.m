function [Qtil,Htil,Ktil,P_all] = UpdateHPIEP (Q,H,K,z,w,pole, poles_old)
% Updating procedure for Hessenberg pencil IEP using [CaMaVaWa20]-criterion
% Given a solution Q,H,K to a HPIEP of size m, the solution to a HPIEP of size m+1
% is constructed, where a new node z, new weight w and new pole p are introduced.
% The updating is done using plane transformations.
%INPUT
%   Q = similarity matrix of HPIEP of size m, such that Q*ZQK=H, for a Z given in HPIEP
%   (H,K) = pencil solving HPIEP of size m
%   z = node to be added HPIEP (= node of inner product)
%   w = weight to be added to HPIEP (= weight of inner product)
%   pole = pole to be added to HPIEP ( = pole of rational Krylov subspace)
%   poles_old = poles of the rational Krylov subspace related to the HPIEP of size m
%     these poles are needed to distinguish the cases that a pole is 0 or infinity
%OUPUT (notation: Ztil = diag(Z,z)
%   Qtil = similarity matrix of size (m+1) such that Qtil*Ztil Qtil Ktil=Htil
%   (Htil,Ktil) = pencil solving HPIEP of size (m+1), with the ratio of their last subdiagonal elements equal to p
% TODO: say where w ends updating
addpath('../../')
m = length(Q);
I = speye(m+1);

%% Embedding
% Embed Q in a larger space, preserving unitarity
Qtil = [Q,zeros(m,1);[zeros(1,m),1]];
% Same for H and K, add a diagonal element
h = z; k = 1;  % must satisfy equality zk = h
Htil = [H,zeros(m,1);[zeros(1,m),h]];
Ktil = [K,zeros(m,1);[zeros(1,m),k]];

%% Orthogonality to weight vector
wold = Q(:,1); % weight vector of HPIEP of size m is given by first column of Q
wtil = [wold; w]; % new weight vector is obtained by adding w to the old one and normalizing

% Construct plane rotation matrix
a = 1/norm(wtil);
b = -w/norm(wtil);

l = 1;

P = I;% P(l,l) = a'; P(m+1,l) = b; P(l,m+1) = -b'; P(m+1,m+1) = a;

P([l,m+1],[l,m+1]) = givens(a,b)';
Qtil = Qtil*P';
Htil = P*Htil;
Ktil = P*Ktil;
% keeps product of all plane rotations working on basis Qtil
P_all = P';
%% Enforcing structure
for l = 2:m
  % The 2x2 block with which we shall work (where element (2,1) will be eliminated)
  Hblock = Htil([l,m+1],[l-1,m+1]);
  Kblock = Ktil([l,m+1],[l-1,m+1]);
  
  temp = Hblock(1,1) * Kblock(2,:) - Kblock(1,1) * Hblock(2,:);
  x1 = temp(1);
  x = temp(2);
  
  %Compute elements of Z = [c', -d';  d, c]
  tau = -x1/x;
  c = sqrt(1/(1+tau*tau'));
  d = c'*tau;
    
  if (Kblock(2,2) == 0 )|| (Kblock(1,1) ~= 0 && (abs(Hblock(2,2))/abs(Kblock(2,2)))>=(abs(Hblock(1,1))/abs(Kblock(1,1))))
    temp = Kblock * [c', -d';  d, c];
  else
    temp = Hblock * [c', -d';  d, c];
  end
  x2 = temp(2,1);
  x = temp(1,1);
  a = sqrt(inv(1+(x2'*x2/(x'*x))));
  b = -x2/x*a;
    

  Pdot = I;% Pdot(l-1,l-1) = c'; Pdot(l-1,m+1) = -d';
       Pdot(m+1,l-1) = d;  Pdot(m+1,m+1) = c;
  Pdot([l-1,m+1],[l-1,m+1]) = givens(c,d)';
  P = I; %P(l,l) = a'; P(l,m+1) = -b'; P(m+1,l) = b; P(m+1,m+1) = a;
  P([l,m+1],[l,m+1]) = givens(a,b)';
  Htil = P*Htil*Pdot;
  Ktil = P*Ktil*Pdot;
  % Set eliminated elements to zero
  Htil(m+1,l-1) = 0;
  Ktil(m+1,l-1) = 0;
  
  Qtil = Qtil*P';
  P_all = P_all*P';
end


%% Introduce pole
mu = Htil(m+1,m);
epsilon = Htil(m+1,m+1);

nu = Ktil(m+1,m);
eta = Ktil(m+1,m+1);
%disp('-------------------------------')
%if mod(m,10)==0
%  m
%  abs(Htil(m:m+1,m:m+1))
%  abs(Ktil(m:m+1,m:m+1))
%end

if isinf(pole)
  rho = -nu/eta;
else
  rho = (mu-nu*pole)/(eta*pole-epsilon);
end
c = sqrt(1/(1+rho'*rho));
d = c'*rho;

Pdot = I; Pdot(m,m) = c'; Pdot(m,m+1) = -d';
       Pdot(m+1,m) = d;  Pdot(m+1,m+1) = c;
Htil = Htil*Pdot;
Ktil = Ktil*Pdot;
if (isinf(pole))
  K(m+1,m) = 0;
end



end
