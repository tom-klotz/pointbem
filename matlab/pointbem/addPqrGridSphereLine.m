function [pqrDataNew,zPoints] = addPqrGridSphereLine(radius, pqrData, spacing)

zPoints = (-radius:spacing:radius)';
numPoints = length(zPoints);
qNew = zeros(numPoints, 1);
xyzNew = [zeros(numPoints,2) zPoints];
R = 0 * qNew;
pqrDataNew = struct('q', [pqrData.q; qNew], 'xyz', [pqrData.xyz; ...
		    xyzNew],'R',[pqrData.R; R]);
