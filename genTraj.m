function genTraj(fname)
  %%% Generate trajectories from the value function saved in 'fname.mat'
  
  %% Settings
  nDisc = 3;
  
  %% Load (
  load(strcat(fname,'.mat'), 'mode', 'X', 'VX', 'Y', 'VY');
  
  %% Process
  [~,dim] = size(Y);
  tMax = Y(1,end);
  xBound = nan(dim-1, 2);
  for i = 1:dim-1
    xBound(i,1) = min(Y(:,i));
    xBound(i,2) = max(Y(:,i));
  end
  scale = (xBound(:,2) - xBound(:,1))';
  
  %% Generate Interp
  if contains(mode, 'kruzkov')
    VX = -log(1-VX);
    VY = -log(1-VY);
  end
  [W, P] = genrbf([X; Y], [VX; VY]);
  
  %% Generate ICs
  nIC = nDisc^dim;
  X0 = nan(nIC, dim);
  prm = permn(linspace(0, 1, nDisc), dim);
  for i = 1:nIC
    for j = 1:dim-1
      X0(i,j) = xBound(j, 1) + prm(i,j)*scale(j);
    end
    j = dim;
    X0(i,j) = prm(i,j)*nDisc*tMax/(nDisc+1);
  end
  
  breakp=1;
  
  %% Integrate trajectories
  
return

%% Functions
  function yD = dydt(t,y)
    
  end

end