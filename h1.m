function c = h1(X)
  %%% Terminal cost - only return non-negative
  s = mySettings();
  x = X(1:end-1);
  t = X(end);
  if isequal(s.problem, 'chilanConway')
    c = 0.5*norm(x)^2/2;
  elseif isequal(s.problem, '1dTargZone')
    c = norm(x) < 2.5;
  elseif isequal(s.problem, 'fancy')
    c = x + 2*sin(4*x);
  elseif isequal(s.problem, 'targZoneTimeCost')
    c = (x < -1) || (x > 1);
  end
return

end