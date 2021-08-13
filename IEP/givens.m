% Computes a complete Givens transformation, like givens.m
% does on my linux computer

function G=givens(x,y)
  [c,s]=cgivens(x,y);
  G=[c,s;-conj(s),c];
