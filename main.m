function main
  %%% Solve for HJB value function approximation via Newton method RBF
  
  
  %% Options
  name = 'timeValTest';
  
  plotAllIters = false;
  % raw_newton, raw_fixpt, kruzkov_newton, OR kruzkov_fixpt
  mode = 'raw_fixpt';
  % special considerations for scenarios with minimum-time-to-capture-type objectives
  s = mySettings();
  objType = s.objType;
  opts = optimoptions('fmincon', 'Display', 'off');
  
  %% Settings
  xBound = [-3, 3];
  tBound = [0, 6];
  uBound = [-1, 1];
  dt = 0.01;

  nSamp = 200;
  
  maxIters = 1000;
  tol = 1e-4;
  
  A = [];
  b = [];
  Aeq = [];
  beq = [];
  nonlcon = [];
  
  
  %% Initialize
  bound = [xBound; tBound];
  scale = (bound(:,2) - bound(:,1))';
  shift = bound(:,1)';
  [dim, ~] = size(bound);

  % Sample terminal surface
  nY = ceil(sqrt(nSamp));
  Y = lhsdesign(nY, dim-1);
  Y(:,dim) = 1;
  for i = 1:nY
    Y(i,:) = scale.*Y(i,:) + shift;
  end
  nEnd = 2^(dim-1);
  prm = dec2bin(nEnd-1:-1:0)-'0' + 1;
  for i = 1:nEnd
    endpoint = nan(dim-1,1);
    for n = 1:dim-1
      endpoint(n,1) = xBound(n, prm(i,n));
    end
    Y(end+1,:) = [endpoint; tBound(2)];
  end
  nY = nY + nEnd;
  if isequal(objType, 'time')
    nTs = ceil(sqrt(nSamp)/2);
    tSamp = linspace(tBound(1), tBound(2), nTs);
    for i = 1:nY
      if h1(Y(i,:)) == 0
        % If y_i is in the target region, its value is known for all t
        yAdd = nan(nTs-1, dim);
        for m = 1:nTs-1
          yAdd(m,:) = [Y(i,1:end-1), tSamp(m)];
        end
        Y = [Y; yAdd];
      end
    end
  end
  [nY, ~] = size(Y);
  % Eval terminal surface
  VY = nan(nY,1);
  for i = 1:nY
    if contains(mode, 'kruzkov')
      VY(i) = 1 - exp(-h1(Y(i,:)));
    else
      VY(i) = h1(Y(i,:));
    end
  end

  % Sample domain
  nX = nSamp - nY;
  X = lhsdesign(nX, dim);
  i = 0;
  while i < nX
    i = i + 1;
    X(i,:) = scale.*X(i,:) + shift;
    if contains(mode, 'fixpt')
      if tBound(2) - X(i,end) < dt
        X(i,end) = tBound(2) - dt;
      end
    end
    if isequal(objType, 'time')
      if h1(X(i,:)) == 0
        % If x is in the target zone, reject and resample
        X(i,:) = rand(1,dim);
        i = i - 1;
      end
    end
  end
  % Initial surface
  for i = 1:nEnd
    endpoint = nan(dim-1,1);
    for n = 1:dim-1
      endpoint(n,1) = xBound(n, prm(i,n));
    end
    X(end+1,:) = [endpoint; tBound(1)];
  end
  nX = nX + nEnd;
  % Init domain = guess  
  if isequal(objType, 'time')
    VX = 0.5*ones(nX, 1);
  else
    for i = 1:nX
      VX(i,:) = interp1(Y(:,1), VY, X(i,1));
    end
  end
  
  avgCtrl = (uBound(:,2) + uBound(:,1))/2;
  
  
  %% Optimize
  % Loop Setup
  targ = 999;
  storeTarg = nan(maxIters, 1);
  iters = 0;
  if contains(mode, 'fixpt')
    nextVX = nan(nX, 1);
  end
  % Until converged
  while targ > tol && iters < maxIters
    iters = iters + 1
    % Solve for RBF weights
    [W, P] = genrbf([X; Y], [VX; VY]);
    if iters == 1 || plotAllIters
      %plotValfun();
    end
    % Calculate HJB constraints and derivatives wrt domain function values
    H = nan(nX,1);
    H_p = nan(nX,nX); 
    for i = 1:nX
      gradRbf = dRbf(X(i,:), W, [X;Y]);
      if isequal(mode, 'kruzkov_newton')
        value = rbf(X(i,:)', W, [X;Y]);
        minfun = @(u) g1(X(i,:), u) + dot(gradRbf(1:end-1), f1(X(i,:), u))...
          - (g(X(i,:), u) - 1)*value;
      elseif isequal(mode, 'kruzkov_fixpt')
        X_delta = @(u) X(i,:)' + dt*[f1(X(i,:), u); 1];
        nonlcon = @(u) boundary_nonlcon(X_delta(u), bound);
        minfun = @(u) rbf(X_delta(u), W, [X;Y])...
          + dt*g1(X(i,:), u)*(1 - VX(i));
      elseif isequal(mode, 'raw_fixpt')
        X_delta = @(u) X(i,:)' + dt*[f1(X(i,:), u); 1];
        nonlcon = @(u) boundary_nonlcon(X_delta(u), bound);
        minfun = @(u) rbf(X_delta(u), W, [X;Y]) + dt*g1(X(i,:), u);
      else  % Raw Newton
        minfun = @(u) g1(X(i,:), u) + dot(gradRbf(1:end-1), f1(X(i,:), u));
      end
      
      [uStar, minVal] = fmincon(minfun, avgCtrl, A, b, Aeq, beq, uBound(:,1), uBound(:,2), nonlcon, opts);
      
      if isequal(mode, 'kruzkov_fixpt')
        nextVX(i) = max(min(minVal, 1), 0);
      elseif isequal(mode, 'raw_fixpt')
        nextVX(i) = minVal;
      end
      
      if contains(mode, 'newton')
        dGradRbfDp = dRbfDy(X(i,:), P, [X;Y]);
        d2vdxdp = dGradRbfDp(1:nX,1:end-1);
      end      
      if isequal(mode, 'kruzkov_newton')
        dValdP_temp = rbfDy(X(i,:), P, [X;Y]);
        dValdP = dValdP_temp(1:nX, :);
        H(i) = value - minVal;
        H_p(i,:) = g(X(i,:), uStar)*dValdP - d2vdxdp*f1(X(i,:), uStar);
      elseif isequal(mode, 'raw_newton')
        dvdt = gradRbf(end,:);
        d2vdtdp = dGradRbfDp(1:nX,end);
        H(i) = dvdt + minVal;
        H_p(i,:) = d2vdtdp + d2vdxdp*f1(X(i,:), uStar);
      end
    end
    % Newton update domain function values
    if contains(mode, 'newton')
      delta = H_p\H;
      VX = VX - delta;
      targ = norm(delta);
    else
      delta = nextVX - VX;
      targ = norm(delta);
      VX = nextVX;
    end
    storeTarg(iters) = targ;
  end
  
  save(strcat(name, '_', mode));
  
  plotValfun();
  breakp=1;

return

%% Functions
  function plotMinfun()
    nP = 100;
    xP = linspace(uBound(1), uBound(2), nP);
    yP = nan(nP, 1);
    for k = 1:nP
      yP(k) = minfun(xP(k));
    end
    figure()
    hold on;
    grid on;
    plot(xP, yP, '-');
    xlabel('Control')
    ylabel('HJB Hamiltonian')
    hold off;
  end

  function plotValfun()
    nP = 50;
    xPlot = linspace(xBound(1), xBound(2), nP);
    tPlot = linspace(tBound(1), tBound(2), nP);
    [XP, TP] = meshgrid(xPlot, tPlot);
    VP = nan(nP, nP);
    for k = 1:nP
      for j = 1:nP
        VP(k,j) = rbf([XP(k,j); TP(k,j)], W, [X;Y]);
      end
    end
    figure()
    hold on;
    grid on;
    mesh(XP, TP, VP);
    plot3([X(:,1); Y(:,1)], [X(:,2); Y(:,2)], [VX; VY], '.');
    xlabel('x')
    ylabel('t')
    zlabel('val')
    hold off;
  end

end