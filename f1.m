function xD = f1(X, u)
  %%% System dynamics - simple velocity control
  x = X(1:end-1);
  t = X(end);
  xD = u;
return

end