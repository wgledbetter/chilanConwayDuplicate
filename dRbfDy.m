function dr = dRbfDy(x, P, X)
  %%% Implementation of derivative of RBF gradient wrt 'training' data Y
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
  
  p = nan(nP, dimX);
  for i = 1:nP
    r_ = x - X(i,:)';
    r = norm(r_);
    dpdr = dPhi(r);
    if r > 0
      drdx = r_/r;
    else
      drdx = 0;
    end
    p(i,:) = dpdr*drdx;
  end
  dr = P'\p;
  
return

end