function A = collocation_mesh(meshData,centroids,normals,areas)

np = size(centroids,1);
% colloc_Laplace and colloc_LPB assume G^L and G^Y are nonnegative.
% this is opposite in sign to the original derivation of Hildebrandt,
% which assumes that G^L is nonpositive, and G^Y is nonnegative.
% HOWEVER Hildebrandt also defines V, K, K' as the negative of the
% singular integrals (see page 102 of his thesis) and D as the regular
% singular integral. so, everything is equal in the end!

[V,K] = colloc_Laplace(meshData,centroids,normals,areas);

A = struct('V',V,'K',K);