function s = mySettings()
  s.phimode = 'thinplatespline';
  s.eps = 0.1;
  % 'chilanConway' OR '1dTargZone' OR 'fancy' OR 'targZoneTimeCost'
  s.problem = 'chilanConway';
  
  if isequal(s.problem, 'targZoneTimeCost')
    s.objType = 'time';
  else
    s.objType = '';
  end
end