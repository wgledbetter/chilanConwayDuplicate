function main2
  %%% Fixed-point iteration for separable two-player games
  
  %% Options
  % raw OR kruzkov
  mode = 'kruzkov';
  opts = optimoptions('fmincon', 'Display', 'off');
  
  %% Settings
  xBound = [0, 4;
            1, 4];
  tBound = [0, 7];
  uBound = [-2*pi, 2*pi];  % Minimizer Control
  vBound = [-1, 1];  % Maximizer Control
  dt = 0.01;
  
  nSamp = 750;
  nEdge = 3;
  
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
%  prm = permn(linspace(0, 1, nEdge), dim-1);
  [nPrm, ~] = size(prm);
  for i = 1:nPrm
    endpoint = nan(dim-1,1);
    for n = 1:dim-1
%      endpoint(n,1) = xBound(n, 1) + prm(i,n)*scale(i,n);
    end
    Y(end+1,:) = [endpoint; tBound(2)];
  end
  nY = nY + nPrm;
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
  for i = 1:nX
    X(i,:) = scale.*X(i,:) + shift;
    if contains(mode, 'fixpt')
      if tBound(2) - X(i,end) < dt
        X(i,end) = tBound(2) - dt;
      end
    end
  end
  % Init domain = guess
  VX = ones(nX, 1);
%  const = mean(VY);
%  VX = const*ones(nX,1);
%  for i = 1:nX
%    VX(i) = interpn(Y(:,1), Y(:,2), VY, X(i,1), X(i,2));
%  end
  
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