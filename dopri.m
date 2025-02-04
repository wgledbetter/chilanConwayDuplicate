function [tOut, yOut, errOut] = dopri(f, tSpan, dt, y0, tol)
  % Dormand-Prince Constant Step-Size Integration
  
  DP = [0,    NaN,        NaN,         NaN,         NaN,      NaN,          NaN,      NaN;
        1/5,  1/5,        NaN,         NaN,         NaN,      NaN,          NaN,      NaN;
        3/10, 3/40,       9/40,        NaN,         NaN,      NaN,          NaN,      NaN;
        4/5,  44/45,      -56/15,      32/9,        NaN,      NaN,          NaN,      NaN;
        8/9,  19372/6561, -25360/2187, 64448/6561, -212/729,  NaN,          NaN,      NaN;
        1,    9017/3168,  -355/33,     46732/5247,  49/176,  -5103/18656,   NaN,      NaN;
        1,    35/384,     0,           500/1113,    125/192, -2187/6784,    11/84,    NaN;
        NaN,  35/384,     0,           500/1113,    125/192, -2187/6784,    11/84,    0;
        NaN,  5179/57600, 0,           7571/16695,  393/640, -92097/339200, 187/2100, 1/40];
      
      
  [tOut, yOut, errOut] = rk(DP, f, tSpan, dt, y0, tol);
end