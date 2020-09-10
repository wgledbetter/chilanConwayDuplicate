function y = rbf(x, W, X)
  %%% Implementation of RBF interpolation given weights W and query x
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
    y = y + W(i,:)*phi(norm(x - X(i,:)'));
  end

return

end