function [collocpts,Zs, areas] = genmeshcolloc(meshData)
% Calculates a list of collocation points by computing panel centroids.
% Panels are stored as 3-D array
% [panel verts, cond num, 0]
% [vert 1 x,y,z]
% [vert 2 x,y,z]
% [vert 3 x,y,z]
% [Other stuff]
%
%
% [vert n x,y,z]
 collocpts(1,3) = 0;
 numpanels = size(meshData.face,1);
 for i=1:numpanels
	numverts = 3;
	panel = [meshData.X(:,i) meshData.Y(:,i) meshData.Z(:,i)];
   [area, collocpt, Z] = calcp(panel);
   collocpts(i,:) = collocpt;
   Zs(i,:) = Z;
   areas(i) = area;
 end

