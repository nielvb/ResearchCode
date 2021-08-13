function [points] = equiPointCircle(k,r,epsilon)
  % Generate the kth point, equidistant on a circle with radius r
  % nodes equidistant gekozen op een cirkel in het complexe vlak
  % Volgorde om de nodes toe te voegen zou telkens het punt moeten zijn met de grootste afstand met de punten die reeds in rekening gebracht zijn.
  
  % k = kth equidistant point on unit circle
  % r = radius of circle
  % epsilon = deviation from its true location, on unit circle
  % 
  % noemer, teller: point = exp(-pi*1i*teller/noemer)
  if nargin<3
    epsilon = 0;
  end
  
  all_k = k;
  points = [];
  
  for i = 1:length(all_k)
    k = all_k(i);
    if k==1
      point = r;
    elseif k==2
      point = -r;
    else
      fk = floor(log2(k-1));
      teller = [1,computeOrder(3:2:2^(fk+1)-1)];
      noemer = 2.^fk;
      p = teller((k-2^fk))/noemer;
      point = r*exp(-pi*1i*(p+epsilon));
    end
    if size(all_k,1)>size(all_k,2)
        points = [points;point];
    else
        points = [points,point];
    end
     
    
  end
