function bem = makeBemYoonLPBMatrices(dielSrf, pqr, epsIn, epsOut, kappaArg)
global kappa;
kappa = kappaArg;

npDiel = length(dielSrf.weights);
nCharges = length(pqr.q);

dielChargeOp = makeSurfaceToChargeOperators(dielSrf, pqr);

Idiel = eye(npDiel);

dielDielLaplace = makeSurfaceToSurfaceLaplaceOperators(dielSrf);
dielDielYukawa =  makeSurfaceToSurfaceYukawaOperators(dielSrf);

A11 = (Idiel/2 + dielDielLaplace.K);
A12 = -dielDielLaplace.V ;
A21 = (Idiel/2 - dielDielYukawa.K);
A22_base = dielDielYukawa.V;

A = [A11 A12;
     A21 (epsIn/epsOut)*A22_base];

B = [dielChargeOp.phiCoul/epsIn;
     zeros(npDiel,nCharges);];

C = 4*pi*[-dielChargeOp.dlpToCharges dielChargeOp.slpToCharges];

bem = struct('B', B, 'A', A, 'C', C,...
	     'dielChargeOp',dielChargeOp,...
	     'dielDielOp',dielDielLaplace,...
	     'dielDielYukawa',dielDielYukawa,...
	     'A11',A11,'A12',A12,...
	     'A21',A21,'A22_base',A22_base);
