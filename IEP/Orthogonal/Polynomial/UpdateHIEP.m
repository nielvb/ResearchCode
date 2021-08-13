function [Qtil,Htil,Ps] = UpdateHIEP(Q,H,z,w)
% Updating procedure for Hessenberg IEP
% Given a solution Q,H to a HIEP of size m, the solution to a HIEP of size m+1
% is constructed, where a new node z, new weight w and new pole p are introduced.
% The updating is done using plane transformations.
%INPUT
%   Q = similarity matrix of HPIEP of size m, such that Q*ZQ=H, for a Z given in HIEP
%   H = Hessenberg matrix solving HIEP of size m
%   z = node to be added to HIEP (= node of inner product)
%   w = weight to be added to HIEP (= weight of inner product)
%OUPUT (notation: Ztil = diag(Z,z)
%   Qtil = similarity matrix of size (m+1) such that Qtil*Ztil Qtil =Htil
%   Htil = Hessenberg matrix solving HIEP of size (m+1)
%   Ps = product of all Ps (plane rotations)

mm = length(Q);
I = speye(mm+1);

%% Embedding
% Embed Q in a larger space, preserving unitarity
Qtil = [Q,zeros(mm,1);[zeros(1,mm),1]];

h = z;
Htil = [H,zeros(mm,1);[zeros(1,mm),h]];

%% Orthogonality to weight vector
wold = Q(:,1); % weight vector of HPIEP of size m is given by first column of Q
wtil_ = [wold; w]; % new weight vector is obtained by adding w to the old one and normalizing

% Construct plane rotation matrix
a = 1/norm(wtil_);
b = -w/norm(wtil_);

l = 1;

P = I;

P([l,mm+1],[l,mm+1]) = givens(a,b)';
Qtil = Qtil*P';
Htil = P*Htil*P';
Ps = P;
%% Enforcing structure
for l = 2:mm
  % The 2x2 block with which we shall work (where element (2,1) will be eliminated)
  Hblock = Htil([l,mm+1],[l-1,mm+1]);
  x2 = Hblock(2,1);
  x = Hblock(1,1);
  a = sqrt(inv(1+(x2'*x2/(x'*x))));
  b = -x2/x*a;
  P = I;% P(l,l) = a'; P(m+1,l) = b; P(l,m+1) = -b'; P(m+1,m+1) = a;

  P([l,mm+1],[l,mm+1]) = givens(a,b)';

  Htil = P*Htil*P';
  % Set eliminated elements to zero
  Htil(mm+1,l-1) = 0;
  
  Qtil = Qtil*P';
  Ps = Ps*P; 
end


end
