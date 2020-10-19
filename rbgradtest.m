function rbgradtest()
  %%% Verify analytical gradients of RBF w.r.t. samples
  % Compare against finite-difference
  
  %% Simple parabola
  if true
    %% Gen RBF
    X = linspace(-1, 1, 11)';
    Y = X.^2;
  
    %% Eval
    x = 0.5;
    
    y = x^2;
    f = theFunction(Y);
    
    fd_g = fd_grad(@theFunction, Y)';
    g = theGradient(Y);
    g - fd_g
    
    fd_gg = fd_grad(@theFunGrad, Y)';
    gg = theFunGradGrad(Y);
    gg - fd_gg
    
    breakp=1;
  end
  
  
return
%% Functions
  function y = theFunction(Y)
    [W1, P1] = genrbf(X,Y);
    y = rbf(x, W1, X);
  end

  function G = theGradient(Y)
    [W1, P1] = genrbf(X,Y);
    G = rbfDy(x, P1, X);
  end

  function yD = theFunGrad(Y)
    [W1, P1] = genrbf(X,Y);
    yD = dRbf(x, W1, X);
  end

  function yDg = theFunGradGrad(Y)
    [W1, P1] = genrbf(X,Y);
    yDg = dRbfDy(x, P1, X);
  end

end