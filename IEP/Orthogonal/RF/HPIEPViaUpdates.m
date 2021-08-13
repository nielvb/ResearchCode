function [Q,H,K] = HPIEPViaUpdates(z,w,p)
% Solves the Hessenbeg pencil IEP using an updating procedure
%
%INPUT:
%     z = nodes
%     w = weights
%     p = poles of rational Krylov subspace (RKS)
%OUPUT:
%     Q = orthogonal basis for RKS with given p
%     (H,K) = Hessenberg pencil such that Q*ZQK = H, with Z = diag(z)

%% Initialize
m = length(z);
Z = diag(z);

% (Trivial) solution to 1x1 IEP
Q = 1;
H = z(1);
K = 1;

%% Solve HPIEP by consecutive updating the solution
normal = norm(w(1));
for i = 2:m
        wtil = [Q(:,1);w(i)/normal];
        [Q,H,K] = UpdateHPIEP_Wa(Q,H,K,z(i),w(i)/normal,p(i-1), p(1:i-2));
        normal = normal*norm(wtil);
end