function [c,s]=cgivens(a,b)
% complex Givens rotation
% determines real c,complex s such that
%   [c       , s][a]=[*]
%   [-conj(s), c][b] [0]

% inspired by blas routine zrotg.f

if a==0,
    c=0;
    s=1;
    return;
end

aa=abs(a);
scale=aa+abs(b);
nrm=scale*sqrt(abs(a/scale)^2+abs(b/scale)^2);
alpha=a/aa;
c=aa/nrm;
s=alpha*conj(b)/nrm;
return;
