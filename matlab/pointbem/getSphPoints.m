function [pntlist, weightlist, angcoords, normals] = getSphPoints(center, ...
																  radius, numpoints);

% method from Park, Bardhan, Makowski, Roux (2009), following
% reference in there

t = (-1+1/numpoints:2/numpoints:1)';
theta = acos(t);
phi = sqrt(pi*numpoints) * asin(t);

normals = [sin(theta).*cos(phi) sin(theta).*sin(phi) cos(theta)];
pntlist = radius * normals;

pntlist = pntlist + ones(numpoints,1) * center;

weightlist = (4 * pi * radius * radius / numpoints) * ...
	 ones(numpoints,1);

angcoords = [theta phi];