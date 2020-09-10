function p = phi(r)
  %%% Radial Basis Function
  s = mySettings();
  eps = s.eps;
  if isequal(s.phimode, 'gaussian')
    mode = @gaussian;
  elseif isequal(s.phimode, 'multiquadric')
    mode = @multiquadric;
  elseif isequal(s.phimode, 'inversequadric')
    mode = @inversequadric;
  elseif isequal(s.phimode, 'inversemultiquadric')
    mode = @inversemultiquadric;
  elseif isequal(s.phimode, 'thinplatespline')
    mode = @thinplatespline;
  elseif isequal(s.phimode, 'bump')
    mode = @bump;
  end
  p = mode(r);
return

%% Functions
  function g = gaussian(r)
    g = exp(-(eps*r)^2);
  end

  function m = multiquadric(r)
    m = sqrt(1 + (eps*r)^2);
  end

  function n = inversequadric(r)
    n = 1/(1 + (eps*r)^2);
  end

  function n = inversemultiquadric(r)
    n = 1/sqrt(1 + (eps*r)^2);
  end

  function t = thinplatespline(r)
    if r == 0
      t = 0;
    else
      t = r^2*log(r);
    end
  end

  function b = bump(r)
    if r < 1/eps
      b = exp(-1/(1 - (eps*r)^2));
    else
      b = 0;
    end
  end

end