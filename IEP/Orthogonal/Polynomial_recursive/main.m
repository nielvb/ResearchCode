close all
clear all
addpath('../../') % to add givens and cgivens to path
n = 10;
m = 30;

z1 = [1:n]';%+rand(n,1)*1i;
z2 = rand(m,1);%+1i;

w1 = ones(n,1);
w2 = rand(m,1);

[Q1,H1] = HIEPViaUpdates(z1,w1);
[Q2,H2] = HIEPViaUpdates(z2,w2);


% checks
norm(Q1(:,1)*norm(w1)-w1)
norm(Q2(:,1)*norm(w2)-w2)

norm(diag(z1)*Q1-Q1*H1)
norm(diag(z2)*Q2-Q2*H2)

% Merge the solutions and create the solution to the joining of both IEPs
[Qtil, Htil]= mergeIEP_HH(Q1,Q2,H1,H2,w1,w2);


% checks
disp('----errors after merge------')
weight_error = norm(Qtil(:,1)*norm([w1;w2])-[w1;w2])
recurrence_error = norm(diag([z1;z2])*Qtil-Qtil*Htil)
struct_error = norm(tril(Htil,-2))
