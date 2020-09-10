function p = dPhi(r)
  %%% Radial Basis Function derivatives
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
    g = -2*r*eps^2*exp(-r^2*eps^2);
  end

  function m = multiquadric(r)
    m = (r*eps^2)/sqrt(r^2*eps^2 + 1);
  end

  function n = inversequadric(r)
    n = -2*r*eps^2/(r^2*eps^2 + 1)^2;
  end

  function n = inversemultiquadric(r)
    n = -r*eps^2/(r^2*eps^2 + 1)^(3/2);
  end

  function t = thinplatespline(r)
    if r == 0
      t = 0;
    else
      t = r + 2*r*log(r);
    end
  end

  function b = bump(r)
    if r < 1/eps
      b = -2*r*x^2*exp(1/(r^2*eps^2 - 1))/(1-r^2*eps^2)^2;
    else
      b = 0;
    end
  end
end