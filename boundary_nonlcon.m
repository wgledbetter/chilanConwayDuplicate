function [c, ceq] = boundary_nonlcon(xDelt, bounds)
  ceq = [];
  c = -[xDelt - bounds(:,1); bounds(:,2) - xDelt];
end