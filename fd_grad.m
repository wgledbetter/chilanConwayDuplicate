function g = fd_grad(func, x)
  eta = 1e-9;
  [r, c] = size(x);
  Idim = max(r,c);
  if c>r
    % Enforce Column Input
    x = x';
  end
  fx = func(x);
  Odim = max(size(fx));
  g = nan(Odim,Idim);
  
  delta = eta*eye(Idim);
  for j = 1:Idim
    xp = x + delta(:,j);
    xn = x - delta(:,j);
    fp = func(xp);
    fn = func(xn);
    dm = (fx-fn)/eta;
    dp = (fp-fx)/eta;
    g(:,j) = (dp+dm)/2;
  end
  
end