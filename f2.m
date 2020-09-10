function xD = f2(X, u, v)
  %%% System dynamics - dolichobrachistrochrone
  x = X(1:end-1);
  t = X(end);
  xD(1,1) = sqrt(x(2))*cos(u) + (v+1)/2;
  xD(2,1) = sqrt(x(2))*sin(u) + (v-1)/2;
end