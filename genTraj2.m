function genTraj2(fname)
  %%% Generate DG trajectories in N dimensions from the value function in 'fname.mat'
  
  %% Settings
  nDisc = 4;
  
  % optimize OR discretize
  trajMode = 'optimize';
  
  intgOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
  
  optOpts = optimoptions('fmincon', 'Display', 'off');
  
  nP = 75;
  
  tStep = 0.01;
  
  %% Load
  load(strcat(fname, '.mat'),...
    'mode', 'X', 'VX', 'Y', 'VY', 'dim', 'xBound', 'tBound', 'uBound', 'vBound',...
    'avgMaxCtrl', 'avgMinCtrl', 'A', 'b', 'Aeq', 'beq', 'dt');
  
  %% Process
  tMax = tBound(2);
  scale = (xBound(:,2) - xBound(:,1))';
  
  %% Generate Interp
  if contains(mode, 'kruzkov')
    VX = -log(1-VX);
    VY = -log(1-VY);
    
    VX(isinf(VX)) = 1000;
    VY(isinf(VY)) = 1000;
  end
  [W, ~] = genrbf([X; Y], [VX; VY]);
  
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
  
  %% Integrate Trajectories
  if isequal(trajMode, 'discretize')
    % Discretization setup
  else
    dydt = @(t,y) dydt_opt(t,y);
  end
  
  paths = cell(nIC, 1);
  parfor i = 1:nIC
    [t, traj, ~] = dopri(dydt, [X0(i,end), tMax], tStep, X0(i,1:end-1), 1e-12);
    paths{i} = [traj, t];
  end
  
  breakp=1;
  
  %% DBC SPECIFIC PLOTTING
  VT = cell(nIC,1);
  for i = 1:nIC
    pth = paths{i};
    [nT, ~] = size(pth);
    val = nan(nT, 1);
    for j = 1:nT
      val(j) = rbf(pth(j,:)', W, [X;Y]);
    end
    VT{i} = val;
  end
  figure()
  ax = axes();
  hold on;
  grid on;
  for i = 1:nIC
    pth = paths{i};
    v = VT{i};
    plot3(pth(:,1), pth(:,2), v, '-', 'LineWidth', 2);
  end
  xlabel('x')
  ylabel('y')
  zlabel('val (t_rem)')
  hold off;
  
  breakp=1;
  
return

%% Functions
  % Dynamics using full optimization
  function yD = dydt_opt(t,y)
    gVal = dRbf([y, t], W, [X;Y]);
    fullfun = @(u,v) -dot(f2([y, t], u, v), -gVal(1:end-1));
    minfun = @(u) fullfun(u, avgMaxCtrl);
    maxfun = @(v) -fullfun(avgMinCtrl, v);
    
    X_delta = @(u, v) y' + dt*f2([y, t], u, v);
    nonlconU = @(u) boundary_nonlcon(X_delta(u, avgMaxCtrl), xBound);
    nonlconV = @(v) boundary_nonlcon(X_delta(avgMinCtrl, v), xBound);
    
    [uStar, ~] = fmincon(minfun, avgMinCtrl, A, b, Aeq, beq, uBound(:,1), uBound(:,2), nonlconU, optOpts);
    [vStar, ~] = fmincon(maxfun, avgMaxCtrl, A, b, Aeq, beq, vBound(:,1), vBound(:,2), nonlconV, optOpts);
    
    yD = f2([y, t], uStar, vStar);
  end

end