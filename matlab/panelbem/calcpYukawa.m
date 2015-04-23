function [phi,phiDouble,phiEdge,phiEdgeDouble] = calcpYukawa(panelOriginal, ...
																  pntsOriginal, ...
																  kappa, panelArea)

% inputs
%  panel: [p1; p2; p3] (column vectors)
%  pnts: [r1; r2; r3 ... r_numpnts] (column vectors)
%  kappa: scalar
% outputs
%  phi: numpnts by 1 (column vector) single layer potential
%  phi_D: double layer potential (not yet implemented)
%  dphidn: gradient of single layer potential (not yet impl)
%  dphi_Ddn: grad of double layer potential (not yet)

numpnts = size(pntsOriginal,2);

% put first vertex at origin
panel = panelOriginal - panelOriginal(:,1) * [1 1 1];
pnts  = pntsOriginal  - panelOriginal(:,1) * ones(1,numpnts);
%for i=1:3
%  fprintf('|pnt - v%d| = %f\n', i, norm(pnts(:,1) - panel(:,i)));
%end

% set up coordinate axes
x = panel(:,2);  x = x/norm(x);
z = cross(panel(:,2), panel(:,3));
z = z/norm(z);
y = cross(z,x);
rotationMat = [x y z];

% put panel in plane
panelPlane = rotationMat' * panel;
pntsPlane  = rotationMat' * pnts;
centroidPlane = mean(panelPlane,2);

% returns a column vector
phi = zeros(numpnts,1);
phiDouble = zeros(numpnts,1);

thresholdMult = 2;
approxPanelRadius = sqrt(panelArea);
%keyboard

% loop
nextIndices = circshift(1:size(panelPlane,2),[1,2]);
for i=1:length(phi)
  panelFinal = panelPlane - [pntsPlane(1:2,i); 0] * ones(1, size(panelPlane,2));
  centroidFinal = centroidPlane - [pntsPlane(1:2,i); 0];
  
  lengthPointToCentroid = norm([0; 0; pntsPlane(3,i)] - centroidFinal);
  if lengthPointToCentroid > thresholdMult * approxPanelRadius
	 expkr = exp(-kappa*lengthPointToCentroid);
	 phi(i) = panelArea * expkr/lengthPointToCentroid;
	 phiDouble(i) = panelArea * (-[0; 0; 1]'*([0; 0; pntsPlane(3,i)]-centroidFinal)) * expkr * (kappa*lengthPointToCentroid+1)/lengthPointToCentroid^3;
  else
	 %  phiEdge = zeros(size(panelPlane,2),1);
	 for j=1:size(panelPlane,2)
		nextJ = nextIndices(j);
		[phiEdge(j,i),phiEdgeDouble(j,i)] = calcpEdgeYukawa(panelFinal(:,j),panelFinal(:,nextJ),...
																  pntsPlane(3,i), kappa);
	 end
	 
	 phi(i) = sum(phiEdge(:,i));
	 phiDouble(i) = sum(phiEdgeDouble(:,i));
  end
end
