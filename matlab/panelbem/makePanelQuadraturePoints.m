function [p, n, w] = makePanelQuadraturePoints(surf, nquad)

% the below code is tailored to interface to
% triangle_unit_set.m. Don't mess with it!
rules = zeros(64,1);
% first column: nquadpoints 
% second column: rule (highest-order rule, if multiple rules with nquadpoints)
mapping = [1    1;
	   3    2;
	   4    4;
	   6    7; % other 6 point rules= precision 3, 7 is precision 4
	   7    9;
	   9   10;
	   12  11;
	   13  12;
	   16  14;
	   64  15;
	   19  17;
	   28  18;
	   37  19];
rules(mapping(:,1))=mapping(:,2);
if rules(nquad) == 0
  fprintf('Error: not a legit nquad!\n');
  p =0; n = 0; w=0;
  return;
end

[xtab,ytab,weight] = triangle_unit_set(rules(nquad));

numpanels = length(surf.areas);
p = zeros(numpanels*nquad,3);
n = 0*p;
w = zeros(numpanels*nquad,1);

for i=1:numpanels
  [R,t,Bg] = makePanelTransform([surf.meshData.X(:,i) ...
		    surf.meshData.Y(:,i) surf.meshData.Z(:,i)]);
  weightTransformed = (surf.areas(i)) * weight; 
  normalTransformed = [0 0 1]*R';
  newPointsTransformed = [[xtab' ytab']*Bg 0*ytab']*R' + repmat(t,[nquad,1]);
  p((i-1)*nquad+1:i*nquad,:) = newPointsTransformed;
  n((i-1)*nquad+1:i*nquad,:) = repmat(normalTransformed,[nquad,1]);
  w((i-1)*nquad+1:i*nquad)   = weightTransformed;
end

% the 1/2 in the weight forces the sum of the weights to = area