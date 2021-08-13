close all
clear

addpath('../../auxiliary')

m =50;


% weights
weights = ones(m,1);

% poles, in complex plane on circle of radius 1.5
rpoles_xi = 3;
poles_xi = equiPointCircle([1:m-1]',rpoles_xi);

rpoles_psi = 3;
poles_psi = equiPointCircle([1:m-1]',rpoles_psi);

% nodes on real line
nodes = equiPointCircle([1:m].',1);
nodes = [real(nodes)+imag(nodes)*0.01*1i];


v = rand(m,1);
w = rand(m,1);
update_proc = true;
if update_proc
  [V,W,T,S,normalv,normalw] = TPIEPViaUpdates(nodes,v,w,poles_xi,poles_psi);
else
  [V,W,T,S,Stil,Ttil] = RationalLanczos_forRatIEP(diag(nodes),v,w,m,poles_xi,ones(size(poles_xi)),poles_psi,ones(size(poles_psi)));
  normalv = v(1)/V(1,1);
  normalw = w(1)/W(1,1);
end


% Check projection
proj_err = norm(diag(nodes)*V*S-V*T)
% Check biorthogonality of basis
orthog_err = norm(W'*V-eye(m))
% Check tridiagonal structure
struct_err = norm(tril(T,-2)+tril(S,-2)+triu(T,2)+triu(S,2)) % Check tridiagonal structure of T and S
 
% Check p
poles_sub_error = poles_xi-diag(T./S,-1);
superDiagRatio = diag(T./S,1);
poles_super_error = conj(poles_psi(1:end-1))-superDiagRatio(2:end);% should be zero
poles_error = max([poles_sub_error;poles_super_error])

% check weights
weigth_error = max([abs(normalv*V(:,1)-v);abs(normalw*W(:,1)-w)])
