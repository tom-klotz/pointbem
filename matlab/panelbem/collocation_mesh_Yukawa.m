function Ay = collocation_mesh_Yukawa(meshData,centroids,normals, ...
												  areas,kappa);

np = size(centroids,1);
[V,K] = colloc_Yukawa(meshData,centroids,normals,areas,kappa);
Ay = struct('V',V,'K',K);
