function operstruct = makeSurfaceToSurfaceOperators(surf)

[V,K] = genPointLaplaceMatrices(surf.points, surf.normals, surf.weights);

if 0
  [V_eigenvecs,V_eigenvals]= eig(V);
  lambda=diag(V_eigenvals);
  [j1,I] = sort(real(lambda),'descend');
  V_eigenvals = lambda(I);
  V_eigenvecs = V_eigenvecs(:,I);
  
  [K_eigenvecs,K_eigenvals]= eig(K);
  lambda=diag(K_eigenvals);
  [j1,I] = sort(real(lambda),'ascend');
  K_eigenvals = lambda(I);
  K_eigenvecs = K_eigenvecs(:,I);
else
  V_eigenvecs = 0; V_eigenvals = 0;
  K_eigenvecs = 0; K_eigenvals = 0;
end

operstruct = struct('V', V, 'K', K, 'V_eigenvals',V_eigenvals, ...
		    'V_eigenvecs',V_eigenvecs, 'K_eigenvals',K_eigenvals, ...
		    'K_eigenvecs',K_eigenvecs);