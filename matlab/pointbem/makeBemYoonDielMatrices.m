function bem = makeBemYoonDielMatrices(dielSrf, pqr, epsIn, epsOut)

npDiel = length(dielSrf.weights);
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargeOperators(dielSrf, pqr);

Idiel = eye(npDiel);

dielDielOp = makeSurfaceToSurfaceLaplaceOperators(dielSrf);

A11 = (Idiel/2 + dielDielOp.K);
A12 = -dielDielOp.V;
A21 = (Idiel/2 - dielDielOp.K);
A22_base = dielDielOp.V;
A = [ A11 A12;
      A21 (epsIn/epsOut)*A22_base;];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);];

C = 4*pi*[-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielOp,...
	     'A11',A11,'A12',A12,...
	     'A21',A21,'A22_base',A22_base);