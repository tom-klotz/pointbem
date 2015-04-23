function [V, K] = colloc_Yukawa(meshData,centroids,normals,areas,kappa)

np = size(meshData.face,1);
for i=1:np
  numverts = 3;
  panel = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
  
  [phi,phiDouble] = calcpYukawa(panel',centroids',kappa, areas(i));
  V(:,i) = phi'/4/pi;
  K(:,i) = phiDouble'/4/pi;
  if mod(i,25)==0
	 fprintf('Done with %d of %d columns.\n',i,np);
  end
end

% the self terms (diagonal entries) always have to be handled by
% themselves.
K = K - diag(diag(K));

