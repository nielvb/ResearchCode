function [Q,H] = solveHIEP(z,w)
% Solves the Hessenbeg IEP using an updating procedure
%
%INPUT:
%     z = nodes
%     w = weights
%OUPUT:
%     Q = orthogonal basis for SKS, with Q e1 = w/norm(w)
%     H = Hessenberg matrix such that Q*ZQ = H, with Z = diag(z)

%% Initialize
m = length(z);
Z = diag(z);

% (Trivial) solution to 1x1 IEP
Q = 1;
H = z(1);

%% Solve HPIEP by consecutive updating the solution
normal = norm(w(1));
for i = 2:m
        wtil = [Q(:,1);w(i)/normal];
        [Q,H] = UpdateHIEP(Q,H,z(i),w(i)/normal);
        normal = normal*norm(wtil);
end

end
