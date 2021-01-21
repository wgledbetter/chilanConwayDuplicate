function c = h2(X)
  %%% Terminal cost - non-negative
  x = X(1:end-1);
  t = X(end);
  %% DBC
  % c = (x(1) ~= 0);  % Incur cost if x1 not zero
  
  %% 1D Simple
  c = (x(1) - x(2))^2;
end