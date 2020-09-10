function [W, P] = genrbf(X, Y)
  %%% Solve for weights of RBF given input X and output Y
  % Input format: Each row a sample
  
  
  %% Setup
  [r,c] = size(X);
  iDim = c;  % Input dimension
  nX = r;
  [r,c] = size(Y);
  oDim = c;  % Output dimension
  nY = r;
  if nX ~= nY
    fail = true
    return
  end
  
  
  %% Build Phi matrix
  P = nan(nX, nX);
  for i = 1:nX
    for j = 1:i
      P(i,j) = phi(norm(X(i,:) - X(j,:)));
      P(j,i) = P(i,j);
    end
  end
  
  
  %% Solve least squares
  W = P\Y;

return

end