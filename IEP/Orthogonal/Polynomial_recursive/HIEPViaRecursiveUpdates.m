function [Q,H] = HIEPViaRecursiveUpdates(z,w)
% Solves Hessenberg IEP with a recursive strategy
% z,w are nodes and weights defining the problem
% Q, H are the unitary and Hessenberg matrix constituting the solution
m = length(z);
% Split up solve and combine, then split up in half as many sets (FFT-style)

Q = eye(m);
H = diag(z);

for i = 1:log2(m)
  for j = 1:2^i:m
    ind = j:j+2^i-1;
    Q1 = Q([ind(1):ind(1)+2^(i-1)-1],[ind(1):ind(1)+2^(i-1)-1]);
    Q2 = Q([ind(1)+2^(i-1):ind(end)],[ind(1)+2^(i-1):ind(end)]);
    
    H1 = H([ind(1):ind(1)+2^(i-1)-1],[ind(1):ind(1)+2^(i-1)-1]);
    H2 = H([ind(1)+2^(i-1):ind(end)],[ind(1)+2^(i-1):ind(end)]);
    
    w1 = w(ind(1):ind(1)+2^(i-1)-1);
    w2 = w(ind(1)+2^(i-1):ind(end));
    
    [Q([ind(1):ind(end)],[ind(1):ind(end)]), H([ind(1):ind(end)],[ind(1):ind(end)])]= mergeIEP_sparse(Q1,Q2,H1,H2,w1,w2);
%    total_counter = total_counter + counter;%%%%%%%%%%%
  end
end



