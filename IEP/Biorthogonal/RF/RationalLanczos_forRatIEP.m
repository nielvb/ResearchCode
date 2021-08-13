function [V,W,T,S,Stil,Ttil] = RationalLanczos_forRatIEP(A,v,w,n,B,L,Lambda,Beta, aux)
% RATLAN computes the biorthogonal bases of rational Krylov subspaces with 
%  given poles and the oblique projections on these subspaces as a pencil.
% Note: also infinite poles are allowed
% Note: small adjustments are made for rational IEP's
%
% This function constructs biorthogonal bases V and W (W*V=I) for two
% rational Krylov subspaces Krat and Lrat. M* denotes the conjugate 
% transpose of some matrix M.
% The subspaces are Krat(A,v,L,B) = span{v,(L(1)A- B(1)I)^{-1}Av,(L(2)A- B(2)I)^{-1}Av}
% and Lrat(A*,w,Beta,Lambda) = span{w,(Beta(1)A*- Lambda(1)I)^{-1}A*v,(Beta(2)A- Lambda(2)I)^{-1}Av}
% In other words B./L are poles of Krat and Lambda./beta of Lrat.
%
% The oblique projection of A onto Krat and orthogonal to Lrat is represented by
% the tridiagonal matrix-pair (T,S) such that W*AVS=T.
% The oblique projection of A* onto Lrat and orthogonal to Krat is represented by
% the tridiagonal matrix-pair (Stil,Ttil) such that V*A*WTtil=Stil.

% INPUT:
%       - A = some (m x m)-matrix 
%       - n = maximum size of subpaces
%       - L, B = variables determining pole i of Krat by their ratio (B(i)/L(i))
%       - Beta, Lambda = variables determining pole i of Lrat by their ratio (Lambda(i)/Beta(i))
%       - v = starting vector for constructing the subspace Krat
%       - w = starting vector for constructing the subspace Lrat
%       - aux = vector containing element (1,2) for T and S 
% Output: with k= min(n+1,m)
%       - V = (m x k)-matrix representing the basis for Krat
%       - W = (m x k)-matrix representing the basis for Lrat
%       - T,S = two tridiagonal matrices of size (k x n) representing the
%       projection of A onto Krat orthogonal to Lrat, i.e., W*AVS=T.
%       - Stil,Ttil = two tridiagonal matrices of size (k x n) representing the
%       projection of A* onto Lrat orthogonal to Krat, i.e., V*A*WTtil=Stil.
  
  
  normal = dot(w,v);
  v = v/sqrt(normal);
  w = w/sqrt(normal)';
  
  V = [v];
  W = [w];
  
  T = zeros(n,1);
  S = zeros(n,1);
  Stil = zeros(n,1);
  Ttil = zeros(n,1);
  
  I = eye(size(A));
  
  G = inv(L(1)*A - B(1)*I);
  d1_c1 = dot(w,G*A*v)/dot(w,G*v);
  vhat = -G*(A-d1_c1*I)*v;
  
  F = inv(Beta(1)*A' - Lambda(1)*I);
  gamma1_delta1 = dot(v,F*A'*w)/dot(v,F*w);
  what = -F*(A'-gamma1_delta1*I)*w;
  
  normal = dot(what,vhat);
  c1 = 1/sqrt(normal);
  delta1 = 1/sqrt(normal)';
  
  vnext = vhat*c1;
  wnext = what*delta1;
  
  T(1,1) = d1_c1*c1;
  T(2,1) = B(1);
  S(1,1) = c1;
  S(2,1) = L(1);
  
  Stil(1,1) = gamma1_delta1*delta1;
  Stil(2,1) = Lambda(1);
  Ttil(1,1) = delta1;
  Ttil(2,1) = Beta(1);
  
  
  V = [V,vnext];
  W = [W,wnext];
  
  %---------------------------------------------------------%
  %---------------------------------------------------------%
  
  
  termination = false;
  if n==length(A)
    n = n-1;
    termination = true;   
  end
  
  
  for i = 2:n
    G = inv(L(i)*A-B(i)*I);
    x = dot(W(:,i),G*A*V(:,i-1));
    y = dot(W(:,i),G*V(:,i-1));
    z = dot(W(:,i),G*A*V(:,i));
    q = dot(W(:,i),G*V(:,i));
    
    xtil = dot(W(:,i-1),G*A*V(:,i-1));
    ytil = dot(W(:,i-1),G*V(:,i-1));
    ztil = dot(W(:,i-1),G*A*V(:,i));
    qtil = dot(W(:,i-1),G*V(:,i));
    
    F = inv(Beta(i)*A'-Lambda(i)*I);
    chi = dot(V(:,i),F*A'*W(:,i-1));
    tau = dot(V(:,i),F*W(:,i-1));;
    eta = dot(V(:,i),F*A'*W(:,i));
    rho = dot(V(:,i),F*W(:,i));
    
    chitil = dot(V(:,i-1),F*A'*W(:,i-1));
    tautil = dot(V(:,i-1),F*W(:,i-1));
    etatil = dot(V(:,i-1),F*A'*W(:,i));
    rhotil = dot(V(:,i-1),F*W(:,i));
    %---------------------------------------------------------%
    if i==2
      betaprev = rand(1);
      lambdaprev = rand(1);  
    else
      betaprev = Beta(i-2);
      lambdaprev = Lambda(i-2);  
    end
    
        
      if betaprev == 0 % psi_{i-2} = infty
        ui = 0;
        di_ai = (y*ztil - ytil*z)/(qtil*z-q*ztil);
        ci_ai = (y+di_ai*q)/z;
        vhat = G*V(:,i-1) - G*(ci_ai*A-di_ai*I)*V(:,i);
      else
        di_ui =  ((lambdaprev/betaprev)'*(y*ztil-ytil*z) + (xtil*z-x*ztil))/(qtil*z-q*ztil);
        ci_ui = ((lambdaprev/betaprev)'*y-x+di_ui*q)/z;
        vhat = -G*(A-(lambdaprev/betaprev)'*I)*V(:,i-1) - G*(ci_ui*A-di_ui*I)*V(:,i);
      end
     
     %---------------------------------------------------------%
    if i==2
      if nargin <9
        lprev = rand(1);
        bprev = rand(1);  
      else
        lprev = aux(1);
        bprev = aux(2);  
      end
    else
      lprev = L(i-2);
      bprev = B(i-2);  
    end
        
      if lprev == 0 % xi_{i-2} = infty
        alphai = 0;
        gammai_mui = (tau*etatil - tautil*eta)/(rhotil*eta-rho*etatil);
        deltai_mui = (tau+gammai_mui*rho)/eta;
        what = F*W(:,i-1) - F*(deltai_mui*A'-gammai_mui*I)*W(:,i);
      else
        gammai_alphai =  ((bprev/lprev)'*(tau*etatil-tautil*eta) + (chitil*eta-chi*etatil))/(rhotil*eta-rho*etatil);
        deltai_alphai = ((bprev/lprev)'*tau-chi+gammai_alphai*rho)/eta;
        what = -F*(A'-(bprev/lprev)'*I)*W(:,i-1) - F*(deltai_alphai*A'-gammai_alphai*I)*W(:,i);
      end
    if(i<length(A))  
     normal = dot(what,vhat);
     vnext = vhat/sqrt(normal);
     wnext = what/sqrt(normal)';
     
     V = [V,vnext];
     W = [W,wnext];
     
     if betaprev == 0 % psi_{i-2} = infty
        ai = 1/sqrt(normal);
        ci = ci_ai*ai;
        di = di_ai*ai;
     else
        ui = 1/sqrt(normal);
        ai = (lambdaprev/betaprev)'*ui;
        ci = ci_ui*ui;
        di = di_ui*ui;
     end
     
     if lprev == 0 % xi_{i-2} = infty   
        mui = 1/sqrt(normal)';
        deltai = deltai_mui*mui;
        gammai = gammai_mui*mui;
     else
        alphai = 1/sqrt(normal)';
        mui = (bprev/lprev)'*alphai;
        deltai = deltai_alphai*alphai;
        gammai = gammai_alphai*alphai;
     end
     
     T(i-1,i) = ai;
     T(i,i) = di;
     T(i+1,i) = B(i);
     S(i-1,i) = ui;
     S(i,i) = ci;
     S(i+1,i) = L(i);
     
     Stil(i-1,i) = mui;
     Stil(i,i) = gammai;
     Stil(i+1,i) = Lambda(i);
     Ttil(i-1,i) = alphai;
     Ttil(i,i) = deltai;
     Ttil(i+1,i) = Beta(i);
   elseif i == length(A)
     if betaprev == 0 % psi_{i-2} = infty
        T(i-1,i) = 1;
        T(i,i) = di_ai;
        S(i-1,i) = 0;
        S(i,i) = ci_ai;
      else
        T(i-1,i) = lambdaprev/betaprev;
        T(i,i) = di_ui;
        S(i-1,i) = 1;
        S(i,i) = ci_ui;
      end
      T = T(1:n,1:n);
      S = S(1:n,1:n);
      
      if lprev == 0 % xi_{i-2} = infty
        Stil(i-1,i) = 1;
        Stil(i,i) = deltai_mui;
        Ttil(i-1,i) = 0;
        Ttil(i,i) = gammai_mui;
      else
        Stil(i-1,i) = bprev/lprev;
        Stil(i,i) = deltai_alphai;
        Ttil(i-1,i) = 1;
        Ttil(i,i) = gammai_alphai;
      end
      Stil = Stil(1:n,1:n);
      Ttil = Ttil(1:n,1:n);
   else
     error('n cannot exceed m')
   end   
  end
  
  if termination && length(A)>2
    m = length(A);
    if Beta(m-1)==0
      am = 0;
      um = 1;
      x = um*V(:,m-1);
      dm = um/(dot(W(:,m-1),A*V(:,m)));
      cm = dm*dot(W(:,m),A*V(:,m));
    else
      psi = Lambda(m-2)/Beta(m-2);
      um = psi';
      x = -(A-psi'*I)*V(:,m-1);
      am = 1;
      dm = am*dot(W(:,m-1),x)/(dot(W(:,m-1),A*V(:,m)));
      cm = dm*dot(W(:,m),A*V(:,m)) - am*dot(W(:,m),x);
    end
    T(m-1,m) = um;
    T(m,m) = cm;
    S(m-1,m) = am;
    S(m,m) = dm;    
  end
     
  