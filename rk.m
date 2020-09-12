function [t, y, err] = rk(bTab, func, tSpan, dt, y0, errMax)
  % CONSTANT STEP SIZE
  % Generic Explicit Runge-Kutta using Butcher Tableau 'bTab'
  calc_error = 0;
  warn_error = 0;
  err = 0;
  
  %% Preliminary
  [rows,cols] = size(bTab);
  stages = cols-1;  % Stages
  if rows>cols  % This means there's an additional row for error management
    calc_error = 1;
    if exist('errMax','var')
      warn_error = 1;
    end
  end
  
  
  dim = max(size(y0));
  
  % Extract 'a', 'b' and 'c' from bTab
  a = bTab(1:cols-1, 2:cols-1);
  b = bTab(cols, 2:cols);
  c = bTab(1:cols-1,1);
  if calc_error
    e = bTab(cols+1, 2:cols);
  end
  
  %% Setup
  t0 = tSpan(1);
  tf = tSpan(2);
  len = tf - t0;
  steps = round(len/dt);
  h = len/steps;  % Just to be sure
  
  t = linspace(t0, tf, steps)';
  y = nan(steps,dim);
  y(1,:) = y0;
  
  if calc_error
    z = nan(steps,dim);
    err = nan(steps,dim);
  end
  %% Integrate
  for i = 1:steps-1
    sum_bk = 0;  % Flush
    k = nan(stages,dim);  % Flush
    for s = 1:stages
      % Calc Ks
      yTerm = 0;  % Flush
      for j = 1:s-1  % Row of Butcher Tableau
        yTerm = yTerm + a(s,j)*k(j,:);
      end
      ys = y(i,:) + h*yTerm;
      ts = t(i) + c(s)*h;
      k(s,:) = func(ts, ys);
      sum_bk = sum_bk + b(s)*k(s,:);
    end
    y(i+1,:) = y(i,:) + h*sum_bk;
    
    if calc_error
      errTerm = 0;
      for s = 1:stages
        errTerm = errTerm + e(s)*k(s,:);
      end
      z(i+1,:) = y(i,:) + errTerm;
      err(i+1,:) = abs( z(i+1,:) - y(i+1,:) );
      if warn_error
        % Check sum of components
        sumErr = sum(err(2:i+1,:));
        violate = sumErr>errMax;
%         if sum(violate)
%           disp('Error Maximum Violated')
%           i
%           errMax
%           sumErr
%         end
      end
    end
    
  end
  
  
  
end