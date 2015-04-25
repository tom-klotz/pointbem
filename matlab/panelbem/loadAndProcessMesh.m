function surfData = loadAndProcessMesh(filename, n1, n2)

areaThreshold = 1e-6;

if nargin < 2
  n1 = 1;
  n2 = 1;
end

meshData = readmesh(filename, n1, n2);

[centroids, normals, areas] = genmeshcolloc(meshData);

I = find(areas > areaThreshold);
areas = areas(I);
centroids = centroids(I,:);
normals = normals(I,:);
meshData.face=meshData.face(I,:);
meshData.X = meshData.X(:,I);
meshData.Y = meshData.Y(:,I);
meshData.Z = meshData.Z(:,I);

surfData = struct('meshData',meshData,'areas',areas, 'centroids',centroids, ...
		  'normals',normals);
