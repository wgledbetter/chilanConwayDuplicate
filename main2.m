function main2
  %%% Fixed-point iteration for separable two-player games
  
  %% Options
  % raw OR kruzkov
  mode = 'raw';
  % time OR nothing
  objType = '';
  opts = optimoptions('fmincon', 'Display', 'off');
  
  %% Settings
  xBound = [-10, 10;
            -10, 10];
  tBound = [0, 10];
  %uBound = [-2*pi, 2*pi];  % Minimizer Control
  uBound = [-1, 1];
  vBound = [-1, 1];  % Maximizer Control
  dt = 0.01;
  
  nSamp = 750;
  nEdge = 5;  % >= 2
  
  maxIters = 1500;
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
  prm = permn(linspace(0, 1, nEdge), dim-1);
  prm(:,dim) = 1;
  [nEnd, ~] = size(prm);
  for i = 1:nEnd
    Y(end+1,:) = scale.*prm(i,:) + shift;
  end
  nY = nY + nEnd;
  if isequal(objType, 'time')
    nTs = ceil(sqrt(nSamp)/2);
    tSamp = linspace(tBound(1), tBound(2), nTs);
    for i = 1:nY
      if h2(Y(i,:)) == 0
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
      VY(i) = 1 - exp(-h2(Y(i,:)));
    else
      VY(i) = h2(Y(i,:));
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
      if h2(X(i,:)) == 0
        % If x is in the target zone, reject and resample
        X(i,:) = rand(1,dim);
        i = i - 1;
      end
    end
  end
  % Initial surface
  prm(:,dim) = 0;
  for i = 1:nEnd
    temp = scale.*prm(i,:) + shift;
    if isequal(objType, 'time')
      if h2(temp) == 0
        % Its already in Y
        continue
      else
        X(end+1,:) = temp;
      end
    end
  end
  [nX, ~] = size(X);
  % Init domain = guess  
  if isequal(objType, 'time')
    VX = 0.5*ones(nX, 1);
  else
%     for i = 1:nX
%       VX(i) = interpn(Y(:,1), Y(:,2), VY, X(i,1), X(i,2));
%     end
  VX = mean(VY)*ones(nX, 1);
  end
  
  avgMinCtrl = (uBound(:,2) + uBound(:,1))/2;
  avgMaxCtrl = (vBound(:,2) + vBound(:,1))/2;
  
  
  %% Optimize
  targ = 999;
  storeTarg = nan(maxIters, 1);
  iters = 0;
  nextVX = nan(nX, 1);
  
  while targ > tol && iters < maxIters
    iters = iters + 1
    
    % Solve for RBF
    [W, P] = genrbf([X; Y], [VX; VY]);
    
    parfor i = 1:nX
      X_delta = @(u, v) X(i,:)' + dt*[f2(X(i,:), u, v); 1];
      nonlconU = @(u) boundary_nonlcon(X_delta(u, avgMaxCtrl), bound);
      nonlconV = @(v) boundary_nonlcon(X_delta(avgMinCtrl, v), bound);
      if isequal(mode, 'kruzkov')
        fullfun = @(u,v) rbf(X_delta(u, v), W, [X;Y])...
          + dt*g2(X(i,:), u, v)*(1 - VX(i));
      elseif isequal(mode, 'raw')
        fullfun = @(u,v) rbf(X_delta(u, v), W, [X;Y])...
          + dt*g2(X(i,:), u, v);
      end
      minfun = @(u) fullfun(u, avgMaxCtrl);
      maxfun = @(v) -fullfun(avgMinCtrl, v);
      
      [uStar, ~] = fmincon(minfun, avgMinCtrl, A, b, Aeq, beq, uBound(:,1), uBound(:,2), nonlconU, opts);
      [vStar, ~] = fmincon(maxfun, avgMaxCtrl, A, b, Aeq, beq, vBound(:,1), vBound(:,2), nonlconV, opts);
      
      breakp=1;

      if isequal(mode, 'kruzkov')
        nextVX(i) = max(min(fullfun(uStar,vStar), 1), 0);
      elseif isequal(mode, 'raw')
        nextVX(i) = fullfun(uStar, vStar);
      end
    end

    delta = nextVX - VX;
    targ = norm(delta);
    VX = nextVX;
    storeTarg(iters) = targ;
  end
  
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
    ylabel('Minimizer HJB Hamiltonian')
    hold off;
  end

  function plotMaxfun()
    nP = 100;
    xP = linspace(vBound(1), vBound(2), nP);
    yP = nan(nP, 1);
    for k = 1:nP
      yP(k) = -maxfun(xP(k));
    end
    figure()
    hold on;
    grid on;
    plot(xP, yP, '-');
    xlabel('Control')
    ylabel('Maximizer HJB Hamiltonian')
    hold off;
  end

end