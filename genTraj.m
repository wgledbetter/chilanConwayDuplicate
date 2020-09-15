function genTraj(fname)
  %%% Generate trajectories from the value function saved in 'fname.mat'
  
  %% Settings
  nDisc = 5;
  nCtrl = 3;
  
  % optimize OR discretize
  trajMode = 'discretize';
  
  intgOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
  
  optOpts = optimoptions('fmincon', 'Display', 'off');
  
  nP = 75;
  
  tStep = 0.01;
  
  %% Load (
  load(strcat(fname,'.mat'),...
    'mode', 'X', 'VX', 'Y', 'VY', 'dim', 'xBound',...
    'uBound', 'avgCtrl', 'A', 'b', 'Aeq', 'beq', 'dt');
  
  %% Process
  tMax = Y(1,end);
  scale = (xBound(:,2) - xBound(:,1))';
  
  %% Generate Interp
  if contains(mode, 'kruzkov')
    VX = -log(1-VX);
    VY = -log(1-VY);
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
  
  breakp=1;
  
  %% Integrate trajectories
  if isequal(trajMode, 'discretize')
    % Create control discretization
    [uDim, ~] = size(uBound);
    prm = permn(linspace(0, 1, nCtrl), uDim);
    nC = nCtrl^uDim;
    uScale = (uBound(:,2) - uBound(:,1))';
    controlOpts = nan(nC, uDim);
    for i = 1:nC
      for j = 1:uDim
        controlOpts(i,j) = uBound(j,1) + prm(i,j)*uScale(j);
      end
    end
    dydt = @(t,y) dydt_disc(t,y);
  else
    dydt = @(t,y) dydt_opt(t,y);
  end
  
  paths = cell(nIC,1);
  parfor i = 1:nIC
%    [t, traj] = ode45(dydt, [X0(i,end), tMax], X0(i,1:end-1));%, intgOpts);
    [t, traj, ~] = dopri(dydt, [X0(i,end), tMax], tStep, X0(i,1:end-1), 1e-12);
    paths{i} = [traj, t];
  end
  
  %% Plot on top of value function
  xPlot = linspace(xBound(1), xBound(2), nP);
  tPlot = linspace(0, tMax, nP);
  [XP, TP] = meshgrid(xPlot, tPlot);
  VP = nan(nP, nP);
  for i = 1:nP
    for j = 1:nP
      VP(i,j) = rbf([XP(i,j); TP(i,j)], W, [X;Y]);
    end
  end
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
  mesh(XP, TP, VP);
  coi = ax.ColorOrderIndex;
  for i = 1:nIC
    ax.ColorOrderIndex = coi;
    pth = paths{i};
    v = VT{i};
    plot3(pth(:,1), pth(:,2), v, '-', 'LineWidth', 2);
  end
  xlabel('x')
  ylabel('t')
  zlabel('val')
  hold off;

  breakp=1;
  
return

%% Functions
  % Dynamics using full optimization
  function yD = dydt_opt(t,y)
    gVal = dRbf([y, t], W, [X;Y]);
    opFun = @(u) -dot(f1([y; t], u), -gVal(1:end-1));
    X_delta = @(u) y' + dt*[f1([y;t], u)];
    nonlcon = @(u) boundary_nonlcon(X_delta(u), xBound);
    
    [uStar, ~] = fmincon(opFun, avgCtrl, A, b, Aeq, beq, uBound(:,1), uBound(:,2), nonlcon, optOpts);
    
    yD = f1([y; t], uStar);
  end

  % Dynamics using discretization selection
  function yD = dydt_disc(t,y)
    gVal = dRbf([y, t], W, [X;Y]);
    opFun = @(u) -dot(f1([y; t], u), -gVal(1:end-1));
    X_delta = @(u) y' + dt*[f1([y;t], u)];
    nonlcon = @(u) boundary_nonlcon(X_delta(u), xBound);
    
    opCase = nan(nC, 1);
    for k = 1:nC
      uk = controlOpts(k,:);
      if any(nonlcon(uk) > 0)
        opCase(k) = Inf;
      else
        opCase(k) = opFun(uk);
      end
    end
    [~, kStar] = min(opCase);
    uStar = controlOpts(kStar,:);
    
    yD = f1([y; t], uStar);
  end

end