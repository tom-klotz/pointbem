function [phi,dphidn] = getPointCoulomb(pqrData, epsIn, points, normals)

numcharges = length(pqrData.q);
numpoints  = size(points,1);
phi = zeros(numpoints,1);
dphidn = phi;

for i=1:numpoints
  for j=1:numcharges
	 if pqrData.q(j) == 0
		continue;
	 end
	 rvec = points(i,:) - pqrData.xyz(j,:);
	 r = norm(rvec);
	 if r < 1e-10
		G = 0; dGdn = 0;
	 else
		G = 1.0/4/pi/r/epsIn;
		dGdn = -dot(rvec,normals(i,:))/4/pi/r^3/epsIn;
	 end
	 phi(i) = phi(i) + pqrData.q(j) * G;
	 dphidn(i) = dphidn(i) + pqrData.q(j) * dGdn;
  end
end
