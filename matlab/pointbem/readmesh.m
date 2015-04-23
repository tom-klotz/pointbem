function vfstruct = readmesh(base,noHeader,noPrune)

vertfile = sprintf('%s.vert', base);
facefile = sprintf('%s.face', base);

if nargin < 2 
  [Xv,Yv,Zv,nx,ny,nz,j1,j2,j3] = textread(vertfile,...
												 '%f %f %f %f %f %f %d %d %d',...
												 'headerlines',3);

  [v1,v2,v3,j1,j2] = textread(facefile,...
										'%d %d %d %d %d',...
										'headerlines',3);
else
  [Xv,Yv,Zv,nx,ny,nz,j1,j2,j3] = textread(vertfile,...
												 '%f %f %f %f %f %f %d %d %d');

  [v1,v2,v3,j1,j2] = textread(facefile,...
										'%d %d %d %d %d');
end

v = [Xv Yv Zv]; 
n = [nx ny nz];
f = [v1 v3 v2];

%%%%%%
if ~ noPrune
thresholdDistance = 1e-4;
thresholdArea     = 1e-6;
faceGood = [];
OKvertices = [];
for i=1:size(f, 1)
  v1 = v(f(i,1),:);
  v2 = v(f(i,2),:);
  v3 = v(f(i,3),:);
  if ((norm(v1-v2) < thresholdDistance) || (norm(v1-v3) < ...
					    thresholdDistance) || (norm(v2-v3) < thresholdDistance))
    faceGood(i) = 0;
  else
    if cross(v2-v1,v3-v1) < thresholdArea
      faceGood(i) = 0;
    else 
      faceGood(i) = 1;
      OKvertices = [OKvertices f(i,1:3)];
    end
  end
end
Iv = unique(OKvertices);
f = f(find(faceGood),:);
%v = v(Iv,:); 
%n = n(Iv,:);
%%%%%
end
for i=1:size(f, 1)
  X(:,i) = v(f(i,1:3),1);
  Y(:,i) = v(f(i,1:3),2);
  Z(:,i) = v(f(i,1:3),3);
end

vfstruct = struct('vert', v, 'face', f, 'normals', n, 'X', X, 'Y', Y, 'Z', Z);
