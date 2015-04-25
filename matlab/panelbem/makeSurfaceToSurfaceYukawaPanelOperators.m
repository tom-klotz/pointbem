function operstruct = makeSurfaceToSurfaceYukawaPanelOperators(surf,kappa,surf2)


if nargin < 3
  [V, K] = colloc_Yukawa(surf.meshData, surf.centroids, surf.normals, ...
			 surf.areas,kappa);
else
  [V, K] = colloc_Yukawa(surf.meshData, surf2.centroids, surf2.normals, ...
			 surf2.areas, kappa);
end

operstruct = struct('V',V, 'K', K);
