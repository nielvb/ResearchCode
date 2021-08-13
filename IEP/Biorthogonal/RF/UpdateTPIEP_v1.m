function [Vtil,Wtil,Ttil,Stil, nuhat,etahat] = UpdateTPIEP_v1 (V,W,T,S,z,v,w,pv,pw, poles_oldv,poles_oldw)
% Updating procedure for tridiagonal pencil IEP
% Given a solution V,W,T,S to a TPIEP of size m, the solution to a TPIEP of size m+1
% is constructed, where a new node z, new weights v,w and new poles pv,pw are introduced.
% The updating is done using eliminators.
%INPUT
%   V,W = biorthonormal matrices of size m, such that W*ZVS=T, for a Z given in TPIEP
%   (T,S) = pencil solving TPIEP of size m
%   z = node to be added TPIEP (= node of bilinear form)
%   v,w = weights to be added to TPIEP (= weights of bilinear form)
%   pv,pw = poles to be added to TPIEP ( = poles of rational Krylov subspace)
%   poles_oldv,poles_oldw = poles of the rational Krylov subspace related to the TPIEP of size m
%     these poles are needed to distinguish the cases that a pole is 0 or infinity
%OUPUT (notation: Ztil = diag(Z,z)
%   Vtil,Wtil = similarity matrix of size (m+1) such that Wtil*Ztil Vtil Stil=Ttil
%   (Ttil,Stil) = pencil solving HPIEP of size (m+1), with the ratio of their last subdiagonal elements equal to p
%   leftprod =  product of all matrices which multiply the embedded (T,S) on the left
% TODO: say where w ends updating

m = length(V);
I = eye(m+1);
%%----------------- Embedding  -----------------%%
% Embed V and W in a larger space, preserving biorthonormality
Vhat = ([V,zeros(m,1);[zeros(1,m),1]]);
What = ([W,zeros(m,1);[zeros(1,m),1]]);
% Same for H and K, add a diagonal element
t = z; s = 1;  % must satisfy equality zk = h
That = ([T,zeros(m,1);[zeros(1,m),t]]);
Shat = ([S,zeros(m,1);[zeros(1,m),s]]);

%% Biorthogonality to weight vector
vold = V(:,1); % weight vector of TPIEP of size m is given by first column of V
wold = W(:,1); % weight vector of TPIEP of size m is given by first column of W
vtil = [vold; v]; % new weight vector is obtained by appending new weight v 
wtil = [wold; w]; 
%%----------------------------------------------%%

% Matrices collecting the products on left and right of the matrices in the pencil
premult = I;
postmult = I;

%%------------- Orthogonality wrt weights  -------------%%
normal_nu = What(:,1)'*vtil;
normal_eta = Vhat(:,1)'*wtil;

% Construct upper and lower matrix which enforce orthogonality L1_dot and R1
%b1 = conj(w);
%a1_dot = -v;%/(normal_eta+v*b1);
%L1_dot = I; L1_dot(m+1,1) = a1_dot/(1+w'*v); L1_dot(1,1) = 1/(1+w'*v);
%R1 = I; R1(1,m+1) = b1;

b1 = conj(w);
a1_dot = -v;
L1_dot = I; L1_dot(m+1,1) = a1_dot;
R1 = I; R1(1,m+1) = b1;
D = I;
D(m+1,m+1)=1/(1-b1*a1_dot);

That = R1*D*L1_dot*That;
Shat = R1*D*L1_dot*Shat;
premult = R1*D*L1_dot*premult;
%%------------------------------------------------------%%
% check weights
%Vhat/premult; ans(:,1)./vtil
%What*premult'; ans(:,1)./wtil

for i = 1:m-1
  if i==1 % this is a special case, deviating from i>1
    % Create colinearity
    c = (That(1,2)*Shat(1,m+1) - Shat(1,2)*That(1,m+1)) / (Shat(1,1)*That(1,m+1) - That(1,1)*Shat(1,m+1));
    C = I; C(1,2)=c;
    That = That*C;
    Shat = Shat*C;
    
    
    a1 = (That(2,1)*Shat(m+1,1) - Shat(2,1)*That(m+1,1)) / (Shat(2,1)*That(m+1,m+1) - That(2,1)*Shat(m+1,m+1));
    L1 = I; L1(m+1,1) = a1;
    That = That*L1;
    Shat = Shat*L1;
    postmult = postmult*C*L1;
%    % Check for colinearity
%    rank([That(2,1),Shat(2,1);That(m+1,1),Shat(m+1,1)],10^-13)
%    rank([That(1,2),That(1,m+1);Shat(1,2),Shat(1,m+1)],10^-13)
    
    % Eliminate the first element in last row and last column of the pencil
    if (1) % TODO: find a suitable criterion
      a2_dot = -That(m+1,1)/That(2,1);
    else
      a2_dot = -Shat(m+1,1)/Shat(2,1);
    end
    if (1)% TODO: find a suitable criterion
      b2_dot = -That(1,m+1)/That(1,2);
    else
      b2_dot = -Shat(1,m+1)/Shat(1,2);
    end
    
    L2_dot = I; L2_dot(m+1,2) = a2_dot;
    R2_dot = I; R2_dot(2,m+1) = b2_dot;
    That = L2_dot*That*R2_dot;
    Shat = L2_dot*Shat*R2_dot;
    % explicitly set eliminates elements to zero ------------
    That(m+1,i) = 0;
    Shat(m+1,i) = 0;
    That(i,m+1) = 0;
    Shat(i,m+1) = 0;
    % ------------------------------------------------------
    premult = L2_dot*premult;
    postmult = postmult*R2_dot;
  else
    Ml = Shat(i+1,i)*[That(i+1,i), That(i+1,m+1); That(m+1,i), That(m+1,m+1)] - That(i+1,i)*[Shat(i+1,i), Shat(i+1,m+1); Shat(m+1,i), Shat(m+1,m+1)];
    Mr = Shat(i,i+1)*[That(i,i+1), That(i,m+1); That(m+1,i+1), That(m+1,m+1)] - That(i,i+1)*[Shat(i,i+1), Shat(i,m+1); Shat(m+1,i+1), Shat(m+1,m+1)];
    % Compute parameters ai and bi of Li and Ri to enforce colinearity
    ai = -Ml(2,1)/Ml(2,2);
    bi = -Mr(1,2)/Mr(2,2);
    
   
   %%%criterion for elimination?
   % Li
%    if (Shat(m+1,m+1) == 0 )|| (Shat(i+1,i) ~= 0 && (abs(That(m+1,m+1))/abs(Shat(m+1,m+1)))>=(abs(That(i+1,i))/abs(Shat(i+1,i))))
%      temp = [Shat(i+1,i), Shat(i+1,m+1); Shat(m+1,i), Shat(m+1,m+1)] * [1,0;ai,1];
%    else
%      temp = [That(i+1,i), That(i+1,m+1); That(m+1,i), That(m+1,m+1)] * [1,0;ai,1];
%    end
   %%%%%%
    
    Li = I; Li(m+1,i) = ai;
    Ri = I; Ri(i,m+1) = bi;
    
    That = Ri*That*Li;
    Shat = Ri*Shat*Li;
    premult = Ri*premult;
    postmult = postmult*Li;
    
    % Eliminate the first element in last row and last column of the pencil
    if (1) % TODO: find a suitable criterion
      anext_dot = -That(m+1,i)/That(i+1,i); 
    else
      ai_dot = -Shat(m+1,i)/Shat(i+1,i);
    end
    if (1)% TODO: find a suitable criterion
      bnext_dot = -That(i,m+1)/That(i,i+1);
    else
      bnext_dot = -Shat(i,m+1)/Shat(i,i+1);
    end
    
    Lnext_dot = I; Lnext_dot(m+1,i+1) = anext_dot;
    Rnext_dot = I; Rnext_dot(i+1,m+1) = bnext_dot;
    That = Lnext_dot*That*Rnext_dot;
    Shat = Lnext_dot*Shat*Rnext_dot;
    if(abs(That(i,i)) <10^-1&& abs(Shat(i,i)) <10^-1 ) || (abs(That(i,i)) >10^1 && abs(Shat(i,i)) >10^1 )
      D = I;
      D(i,i) = 10^(-log10(That(i,i)));
      That = That*D;
      Shat = Shat*D;
    end
    
    % explicitly set eliminates elements to zero ------------
    That(m+1,i) = 0;
    Shat(m+1,i) = 0;
    That(i,m+1) = 0;
    Shat(i,m+1) = 0;
    % ------------------------------------------------------
    
    premult = Lnext_dot*premult;
    postmult = postmult*Rnext_dot;
    
    if abs(That(m+1,m+1)<10^-8) || abs(Shat(m+1,m+1)<10^-8)
      disp('###### Warning: close to numerical breakdown ########################')
    end
    
  end
end

%--------------------Introduce poles------------------------%
xi = pv;
if isinf(xi)
  am = - Shat(m+1,m) / Shat(m+1,m+1);
else
  % split pole in fraction, currently: mu = pole and nu =1, pv = mu/nu
  nu = 1;
  mu = xi*nu;
  am = (mu*Shat(m+1,m) - nu*That(m+1,m)) / (nu*That(m+1,m+1) - mu*Shat(m+1,m+1));
end
Lm = I; Lm(m+1,m) = am;

if m>1
  psi = poles_oldw(end);
  if isinf(psi)
    bm = - Shat(m,m+1) / Shat(m+1,m+1);
  else
    % split pole in fraction, currently: alpha = pole and beta =1, pw = alpha/beta
    beta = 1;
    alpha = psi*beta;
    bm = (alpha'*Shat(m,m+1) - beta'*That(m,m+1)) / (beta'*That(m+1,m+1) - alpha'*Shat(m+1,m+1));
  end
    Rm = I; Rm(m,m+1) = bm;
else 
   Rm = I;
end

%-----------------------------------------------------------%
Ttil = Rm*That*Lm;
Stil = Rm*Shat*Lm;

premult = Rm*premult;
postmult = postmult*Lm;

Vtil = Vhat/premult;
Wtil = What*premult';


nuhat = Wtil(:,1)'*vtil;
etahat = Vtil(:,1)'*wtil;
end
