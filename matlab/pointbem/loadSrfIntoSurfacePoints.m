function surfdata = loadSrfIntoSurfacePoints(srf)

[meshBase,rootDir] = readsrf(srf);
meshData = readmesh(meshBase,1,0);
[panelCentroids,panelNormals,panelAreas] = genmeshcolloc(meshData);

vertexWeights = zeros(size(meshData.vert,1),1);
for i=1:size(meshData.face,1)
  vertexWeights(meshData.face(i,1:3)) = vertexWeights(meshData.face(i,1:3)) + panelAreas(i)/3.0;
end
thresholdDistance = 1e-6;
thresholdWeight = 1e-4;
eliminateList = ones(length(vertexWeights),1);
keyboard
for i=1:length(vertexWeights)
  if vertexWeights(i) < thresholdWeight
    eliminateList(i) = 0;
  else
    for j=i+1:length(vertexWeights)
      if (norm(meshData.vert(i,:)-meshData.vert(j,:)) < ...
	  thresholdDistance)
	eliminateList(j) = 0;
	vertexWeights(i) = vertexWeights(i) + vertexWeights(j);
	vertexWeights(j) = 0;
      end
    end
  end
end
meshData.vert = meshData.vert(find(eliminateList),:);
vertexWeights = vertexWeights(find(eliminateList));
meshData.normals = meshData.normals(find(eliminateList),:);
surfdata = struct('points',meshData.vert,'weights',vertexWeights, ...
		  'anglecoords', 0,'normals',meshData.normals);

% to plot: 
% v = meshData.vert;
% n = meshData.normals;
% quiver3(v(:,1),v(:,2),v(:,3),n(:,1),n(:,2),n(:,3));
