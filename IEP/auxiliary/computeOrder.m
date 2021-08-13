 function order = computeOrder(p)
   
   if length(p) == 1
     order = p;
   else
     ind = (length(p)+1)/2;
     order = [p(ind)];
     orderleft = computeOrder(p(1:ind-1));
     orderright = computeOrder(p(ind+1:end));
     c = 2^0;
     while c < ind
       subsetLeft = orderleft(c:2*c-1);
       subsetRight = orderright(c:2*c-1);
       % merge subsets correctly
       orderedSubset(1:2:length(subsetLeft)*2) = subsetLeft;
       orderedSubset(2:2:length(subsetLeft)*2) = subsetRight;
       order = [order, orderedSubset];
       c = 2*c;
     end
     
   end
   