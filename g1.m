function c = g1(X, u)
  %%% Integral cost - only return non-negative
  s = mySettings();
  x = X(1:end-1);
  t = X(end);
  if isequal(s.problem, 'chilanConway')
    c = norm(u)^2/2;
  elseif isequal(s.problem, '1dTargZone')
    c = 0;
  end
return

end