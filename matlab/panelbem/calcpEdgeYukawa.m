function [phiEdge,phiEdgeDouble] = calcpEdgeYukawa(p1, p2, z, kappa)

phiEdge = 0;
phiEdgeDouble = 0;

p21 = p2 - p1;
lengthp21 = norm(p21);
p21normalized = p21/lengthp21;
orthog = cross([0 0 1]',p21normalized);

alpha = -(p21'*p1)/lengthp21^2;
rOrthog = p1 + alpha * p21;
distanceToEdge = norm(rOrthog);
whichSideVec = cross(p2-p1,-p1);

rotateToVertLine = [orthog'; p21normalized'; 0 0 1];
p1new = rotateToVertLine * p1;
if p1new(1) < 0
  p21normalized = -p21normalized;
  orthog = -orthog;
  rotateToVertLine = [orthog'; p21normalized'; 0 0 1];
  p1new = rotateToVertLine * p1;
end
p2new = rotateToVertLine * p2;
rOrthognew = rotateToVertLine * rOrthog;
x = p1new(1);

if (((p1new(2) > 0) && (p2new(2) < 0)) || ...
	 ((p1new(2) < 0) && (p2new(2) > 0)) )

  [phiEdge1,phiEdgeDouble1] = logIntDoLine(kappa, z, x, 0,p1new(2),1);
  [phiEdge2,phiEdgeDouble2] = logIntDoLine(kappa, z, x, p2new(2),0,1);

  phiEdge = phiEdge1+phiEdge2;
  phiEdgeDouble = phiEdgeDouble1 +phiEdgeDouble2;

else

  [phiEdge,phiEdgeDouble] = logIntDoLine(kappa, z, x, p1new(2),p2new(2));
  phiEdge = -phiEdge;
  phiEdgeDouble = -phiEdgeDouble;
end
