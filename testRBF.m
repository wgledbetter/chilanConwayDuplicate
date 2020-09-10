function testRBF
  %%% Test my implementation of RBF
  
  %% Settings
  nSamp = 50;
  nPlot = 100;
  xPlot = linspace(-2.5, 2.5, nPlot);
  mu = 0;
  sig = 1;
  
  %% 1D
%   dimX = 1;
%   X = rand(nSamp, dimX);
% %  X = lhsnorm(mu, sig, nSamp);
%   Y = seedFun(X);
%   
%   % Solve
%   [W, P] = genrbf(X, Y);
%   
%   % Predict and Plot
%   for i = 1:nPlot
%     yPlot(i,1) = rbf(xPlot(i), W, X);
%   end
%   
%   figure()
%   hold on;
%   grid on;
%   plot(X, Y, '.', 'MarkerSize', 10);
%   plot(xPlot, yPlot, '-');
%   hold off;
%   
%   
%   breakp=1;
  
  
  %% 2D
%   dimX = 2;
%   X = rand(nSamp, dimX);
%   X = lhsnorm([mu mu], diag([sig sig]), nSamp);
%   Y = seedFun(X);
%   
%   % Solve
%   [W, P] = genrbf(X, Y);
%   
%   % Predict and Plot
%   [xMesh1, xMesh2] = meshgrid(xPlot, xPlot);
%   for i = 1:nPlot
%     for j = 1:nPlot
%       yMesh(i,j) = rbf([xMesh1(i,j), xMesh2(i,j)], W, X);
%     end
%   end
%   
%   figure()
%   hold on;
%   for i = 1:nSamp
%     plot3(X(i,1), X(i,2), Y(i), '.', 'MarkerSize', 10)
%   end
%   mesh(xMesh1, xMesh2, yMesh);
%   hold off;
%   
%   
%   breakp=1;
  

  %% 1D Gradients
  dimX = 1;
  X = rand(nSamp, dimX);
  Y = seedFun(X);
  
  % Solve
  [W, P] = genrbf(X, Y);
  
  % Calculate gradients at random point
  x = rand(1,dimX);
  gradX = dRbf(x, W, X);
  gradXdY = dRbfDy(x, P, X);
  
  breakp=1;
  
  
  %% 3D Gradients
  dimX = 3;
  X = rand(nSamp, dimX);
  Y = seedFun(X);
  
  % Solve
  [W, P] = genrbf(X, Y);
  
  % Calculate gradients at random point
  x = rand(1,dimX);
  gradX = dRbf(x, W, X);
  gradXdY = dRbfDy(x, P, X);
  
  breakp=1;

return

  function y = seedFun(x)
    [nX,~] = size(x);
    for i = 1:nX
%      y(i,1) = norm(x(i,:))^2;
      y(i,1) = sin(norm(3*x(i,:)));
    end
  end

end