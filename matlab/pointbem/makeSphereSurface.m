function surfdata = makeSphereSurface(origin, radius, numpoints)

[pnts,wts,angs,normals] = getSphPoints(origin, radius, numpoints);

surfdata = struct('points',pnts,'weights',wts, 'anglecoords', ...
		  angs,'normals',normals);