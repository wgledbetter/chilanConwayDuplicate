function dr = rbfDy(x, P, X)
  %%% Implementation of RBF derivative wrt 'training' data Y
  % This could be implemented without passing in P, but its faster this way
  % and removes redundancy
  [r,c] = size(x);
  if c > r
    x = x';
  end
  dimX = max([c r]);
  [nW,~] = size(P');
  [r,c] = size(X);
  nP = r;
  dimData = c;
  if (nP ~= nW) || (dimX ~= dimData)
    fail = true
    return;
  end
  
  p = nan(nP,1);
  for i = 1:nP
    p(i) = phi(norm(x - X(i,:)'));
  end
  dr = P'\p;
  
return

end