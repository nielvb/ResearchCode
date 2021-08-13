function [V,W,Tv,Sv, normalv,normalw] = TPIEPViaUpdates(z,v,w,pv,pw)
% Solves the tridiagonal pencil IEP using an updating procedure
%
%INPUT:
%     z = nodes
%     v,w = weights
%     pv,pw = poles of rational Krylov subspaces (RKS)
%OUPUT:
%     Q = orthogonal basis for RKS with given p
%     (H,K) = Hessenberg pencil such that Q*ZQK = H, with Z = diag(z)

%% Initialize
m = length(z);
Z = diag(z);

% (Trivial) solution to 1x1 IEP
normal = dot(w(1),v(1));
normalv = sqrt(normal);
normalw = normalv';
V = v(1)/normalv;
W = w(1)/normalw;
Tv = z(1);
Sv = 1;

%Vf = 1;
%Wf = 1;
%Tvf = z(1);
%Svf = 1;

%% Solve HPIEP by consecutive updating the solution
for i = 2:m
        [V,W,Tv,Sv,nuhat,etahat] = UpdateTPIEP_v1(V,W,Tv,Sv,z(i),v(i)/normalv,w(i)/normalw,pv(i-1),pw(i-1), pv(1:i-2),pw(1:i-2));
%        [Vf,Wf,Tvf,Svf,~] = UpdateTPIEP (Vf,Wf,Tvf,Svf,z(i),v(i),w(i),pv(i-1),pw(i-1), pv(1:i-2),pw(1:i-2));
        normalv = normalv*nuhat;
        normalw = normalw*etahat;
        if abs(normalv)<10^-8 || abs(normalw)<10^-8
          disp('***********bilinear form vanishing**********************')
        end
end