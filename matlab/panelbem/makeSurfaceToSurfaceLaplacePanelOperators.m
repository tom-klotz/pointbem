function operstruct = makeSurfaceToSurfaceLaplacePanelOperators(surf, ...
						  surf2)

if nargin < 2
  [V, K, Kp] = colloc_Laplace(surf.meshData, surf.centroids, surf.normals, ...
			      surf.areas);
else
  [V, K, Kp] = colloc_Laplace(surf.meshData, surf2.centroids, ...
			      surf2.normals, surf2.areas);
end

operstruct = struct('V',V, 'K', K, 'Kp', Kp);