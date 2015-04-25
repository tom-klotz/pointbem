function [G,dGdn] = coulombMatrixCollocation(meshData,...
															centroids,normals, ...
															chargeLocations)
np = size(centroids,1);
ncharges = size(chargeLocations,1);

G = zeros(np,ncharges);
dGdn = G;

% coded badly. apologies.
for i=1:np
  for j=1:ncharges
	 r = centroids(i,:) - chargeLocations(j,:);
	 normr = norm(r);
	 G(i,j) = 1/normr; % see eq 5.103 on Hildebrandt thesis pg 107
	 dGdn(i,j) = -r*normals(i,:)'/(normr^3);
  end
end

G = G/4/pi;
dGdn = dGdn/4/pi;