function xD = f2(X, u, v)
  %%% System dynamics
  x = X(1:end-1);
  t = X(end);
  %% Dolichobrachistochrone
  % xD(1,1) = sqrt(x(2))*cos(u) + (v+1)/2;
  % xD(2,1) = sqrt(x(2))*sin(u) + (v-1)/2;
  
  %% 1D Simple Motion Pursuit
  xD(1,1) = u(1);  % Pursuer
  xD(2,1) = v(1);  % Evader
  
end