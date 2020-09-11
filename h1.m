function c = h1(X)
  %%% Terminal cost - only return non-negative
  s = mySettings();
  x = X(1:end-1);
  t = X(end);
  if isequal(s.problem, 'chilanConway')
    c = 0.5*norm(x)^2/2;
  elseif isequal(s.problem, '1dTargZone')
    c = norm(x) < 2.5;
  end
return

end