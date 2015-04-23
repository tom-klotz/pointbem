function pqrData = makeSphereChargeDistribution(radius, numCharges, ...
						gridDistance, minDistance)

if nargin < 4
  minDistance = gridDistance;
end

maxChargeValue = 0.85;

mygrid = -radius:gridDistance:radius;
[x,y,z]=meshgrid(mygrid,mygrid,mygrid);

xyz = [reshape(x,numel(x),1) reshape(y,numel(y),1) reshape(z, ...
																  numel(z),1)];
gridOk = zeros(numel(x),1);
for i=1:size(xyz,1)
  if norm(xyz(i,:)) < radius - minDistance
	 gridOk(i) = 1;
  end
end

xyz = xyz(find(gridOk>0),:);
q   = 2*maxChargeValue*(rand(size(xyz,1),1) - 0.5);

if ((numCharges < 0) || (numCharges > length(q)))
  fprintf('Info: using all the grid charges in makeSphereChargeDistribution\n');
else 
  I = randperm(length(q)); 
  I = I(1:min(length(I),numCharges));
  q = q(I); 
  xyz = xyz(I,:);
end

pqrData = struct('q',q, 'xyz', xyz,'R', 0*q);