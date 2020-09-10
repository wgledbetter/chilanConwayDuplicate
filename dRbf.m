function y = dRbf(x, W, X)
  %%% Gradient of RBF interpolation given weights W and query x
  [r,c] = size(x);
  if c > r
    x = x';
  end
  dimX = max([c r]);
  [nW,~] = size(W);
  [r,c] = size(X);
  nP = r;
  dimData = c;
  if (nP ~= nW) || (dimX ~= dimData)
    fail = true
    return;
  end
  
  y = 0;
  for i = 1:nP
    r_ = x - X(i,:)';
    r = norm(r_);
    dpdr = dPhi(r);
    if r > 0
      drdx = r_/r;
    else
      drdx = 0;
    end
    y = y + W(i,:)*dpdr*drdx;
  end
return

end