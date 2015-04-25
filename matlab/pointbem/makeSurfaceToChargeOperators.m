function chargesurfop = makeSurfaceToChargeOperators(surf, pqr)
nq = length(pqr.q);
np = length(surf.weights);

phiCoul = zeros(np,nq); 
dphidnCoul = 0 * phiCoul;
slpToCharges = 0 * phiCoul';
dlpToCharges = 0 * slpToCharges;

for i=1:np
  for j=1:nq
    rvec = surf.points(i,:) - pqr.xyz(j,:);
    r = norm(rvec);

    if r < 1e-10
      G = 0; dGdn = 0;
    else
      G = 1.0/4/pi/r;
      dGdn = -dot(rvec,surf.normals(i,:))/4/pi/r^3;
    end
    phiCoul(i,j) = G;
    dphidnCoul(i,j) = dGdn;
    slpToCharges(j,i) = G * surf.weights(i);
    dlpToCharges(j,i) = dGdn * surf.weights(i);
  end
end

chargesurfop = struct('phiCoul', phiCoul,...
		      'dphidnCoul',dphidnCoul,...
		      'slpToCharges',slpToCharges,...
		      'dlpToCharges',dlpToCharges);


